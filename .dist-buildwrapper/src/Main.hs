{-# LANGUAGE ScopedTypeVariables, RankNTypes #-}
module Main where

import qualified Data.Map.Strict as Map 

import SolutionDomain
import CalculusUtils
import Data.List
import Control.Monad.Reader as Reader

--{-
import Control.Monad.Par as Par
import qualified Graphics.Gnuplot.Graph as Graph
import qualified Graphics.Gnuplot.Frame as Frame
import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import qualified Graphics.Gnuplot.Frame.Option as Opt
import qualified Graphics.Gnuplot.Advanced as GP
import qualified Graphics.Gnuplot.Plot.TwoDimensional as Plot2D
import qualified Graphics.Gnuplot.Graph.TwoDimensional as Graph2D
import Graphics.Gnuplot.Terminal.PNG as PNG 

instance Par.NFData Property 
instance Par.NFData Position

 
defltOpts :: Graph.C graph => Opts.T graph
defltOpts = Opts.key False Opts.deflt
    
plotDomain grid scale= 
    let xsize = length (head grid)-1 
        ysize = length  grid  -1
        printSize = show (ysize * 5) ++ "," ++ show (xsize * 5)
    in Frame.cons 
            (   Opts.add ( Opt.custom "terminal pngcairo size" printSize) [printSize] $ 
                scale $
                Opts.sizeRatio (fromIntegral xsize / fromIntegral ysize)  defltOpts 
            ) $
        Plot2D.function Graph2D.image 
        (liftM2 (,)  [0..ysize] [0..xsize]) 
            $ \(x,y) -> (grid!!x)!!y
             
plotDomainLog grid = plotDomain grid $ Opts.add ( Opt.custom "logscale cb"  "") []

plotDomainLinear grid = plotDomain grid id
---}

heatConductivityK:: Double
heatConductivityK = 0.025

orthogonalSides:: Side ->[Side]
orthogonalSides side = case side of
    East -> [North, South, Top, Bottom]
    West -> [North, South, Top, Bottom]
    North -> [East, West, Top,Bottom]
    South -> [East, West, Top,Bottom]
    Top -> [East,West,North,South]
    Bottom -> [East,West,North,South]
    Now -> [Now,Prev]
    Prev -> [Now,Prev]
    Center -> [Center]

addTerm:: [Term a]->Term a->[Term a]
addTerm terms term = terms ++ [term]

addTerms:: [Term a]-> [Term a]-> [Term a]
addTerms terms1 terms2 = case terms2 of
    ([x]) -> addTerm terms1 x
    (x:xs) -> addTerms (addTerm terms1 x) xs
    _ -> terms1
        
testEquation:: (Num a, Fractional a ) => Equation (Term a )
testEquation = 
    Equation 
        --[Unknown 2, Constant 4, Constant (-0.232)]
        [Unknown 2, SubExpression (Expression [Constant 2, Constant 2, Unknown 0.025]), Constant (-0.232)] 
        ( addTerms [Unknown (-0.025),Unknown (-0.05)] [Constant 2, Constant 3 ] )
        U

--continuity::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
--continuity::Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = do
    env <- ask
    return $ let integrate = integUnknown env Density 
        in Equation
            ((integrate Temporal [runReader drho_dt env] >>= integrate Spatial) 
              : integrateTerms integrate env ( divergenceWithProps [Density]) )
            [ const [Constant 0]] 
            Density
    
--uMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
--uMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
uMomentum = do
    env <- ask
    return $ let integrate = integUnknown env U
        in Equation
            (  [ integrate Temporal [runReader drhou_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density, U])
                 ++ [integrate Spatial [runReader dp_dx env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [  divGrad [Mew,U] 1
                  ,divergence [ dmewu_dx, dmewv_dx, dmeww_dx ]
                  ,map (\d-> ddf_ d (-2/3) X) (divergenceWithProps [Mew]) 
                ] 
            )
            U

--vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
vMomentum = do
    env <- ask
    return $ let integrate = integUnknown env V
        in Equation
            ([ integrate Temporal [runReader drhov_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density, V])
                ++ [integrate Spatial [runReader dp_dy env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,V] 1 
                  ,divergence [ dmewu_dy, dmewv_dy, dmeww_dy ] 
                  ,map (\d-> ddf_ d (-2/3) Y) (divergenceWithProps [Mew]) 
                ]
            )
            V   

--wMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
wMomentum =  do
    env <- ask
    return $ let integrate = integUnknown env W 
        in Equation
            ([ integrate Temporal [runReader drhow_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density , W] ) 
                ++[integrate Spatial [runReader dp_dz env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,W] 1 
                  ,divergence [ dmewu_dz, dmewv_dz, dmeww_dz ] 
                  ,map (\d-> ddf_ d (-2/3) Z) (divergenceWithProps [Mew]) 
                ]
            )
            W      

dye = do
    env <- ask
    return $ let integrate = integUnknown env Dye 
        in Equation
            ([ integrate Temporal [runReader drhodye_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithProps [Density , Dye] ) 
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,Dye] 0.001  ]
            )
            Dye

--energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
energy =  do
    env <- ask
    return $ let integrate = integUnknown env Temperature
        in Equation
            ([ integrate Temporal [runReader drhoT_dt env] >>= integrate Spatial ]
                ++ integrateTerms integrate env (divergenceWithPropsFactor [Density,Temperature] specificHeatCv ) 
                ++ integrateTerms integrate env (divergenceWithProps [Pressure] )
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Temperature] heatConductivityK 
                  , squareDerivative [Mew,U] (8/3) X 
                  , squareDerivative [Mew,V] (8/3) Y   
                  , squareDerivative [Mew,W] (8/3) Z
                  , squareDerivative [Mew,U] 1 Y
                  , squareDerivative [Mew,V] 1 X
                  , squareDerivative [Mew,U] 1 Z
                  , squareDerivative [Mew,W] 1 X
                  , squareDerivative [Mew,V] 1 Z
                  , squareDerivative [Mew,W] 1 Y
                  , pairedMultipliedDerivatives [Mew,U][V] X Y
                  , pairedMultipliedDerivatives [Mew,U][W] X Z
                  , pairedMultipliedDerivatives [Mew,W][V] Z Y                         
                ]
            )
            Temperature   
    
getDiscEqInstance :: forall a t.Equation (t -> [a]) -> t -> Equation a
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) $! l) (concatMap (\t -> t pos) $! r) up

applyDiffEq (eq, saveAtNextTime,getPos) env =
    runPar $ parMap
    --map
        (\pos -> 
            let discEquation = getDiscEqInstance (runReader eq env) pos 
                solvedProperty = unknownProperty discEquation
                newValue = solveUnknown discEquation pos env
            in ( pos, solvedProperty, newValue, saveAtNextTime)
        ) (getPos env)
  

applyResults :: [(Position, Property, Double, Bool)] -> Bool -> ValSet Double -> ValSet Double
applyResults res pushTime vs = 
    let (ValSet p v av sl tl) = foldl' 
            (\prevEnv (pos,property,newVal,saveAtNextTime)  ->
                let newPos = if saveAtNextTime 
                        then advanceTime pos 
                        else pos
                    newValAdded = setVal prevEnv newPos property newVal -- really the intermediate value has been added, wil be corrected next time step
                in if timeLevelAbsolute newValAdded > 0 
                    then 
                        let prevVal = prop Nondirectional property (offsetPosition pos Prev) Center newValAdded
                            correctedCurrVal = prevVal + 0.5 * ( newVal - prevVal) 
                        in setVal newValAdded pos property correctedCurrVal
                    else newValAdded -- there are not enough time steps yet to do Heun's Method
            ) vs res
    in ValSet 
        (if pushTime then map advanceTime p else p) -- this is what determines the list of positions (inclduing time levels) of the "current" level
        v av sl (if pushTime then tl+1 else tl)

calcSteps = [ 
    (continuity, True, allPositionsCurrentTime )
    ,(uMomentum, True, calculatedPositions )
    ,(vMomentum,  True, calculatedPositions )
    ,(wMomentum,  True, calculatedPositions )
    ,(energy,  True, allPositionsCurrentTime )  
    ,(dye, True, calculatedPositions )
    ] 

supportCalcSteps = []

allPositionsCurrentTime :: forall a. ValSet a -> [Position]
allPositionsCurrentTime env = 
    let curTimeLevel = timePos $ head $ calculatedPositions env
    in map (\x-> modifyPositionComponent x Time curTimeLevel ) makeAllPositions 

runSingleStep prev rerun = 
    let runSteps steps vs= concat $!
            --runPar $ parMap
            map (\step -> snd $ foldl' 
                    (\(prevEnv,prevResList) _ ->
                        let nextResList = applyDiffEq step nextEnv
                            nextEnv = applyResults prevResList False prevEnv 
                        in (nextEnv, nextResList)
                    ) (vs , applyDiffEq step vs) [1..2]  
                ) steps  
    in foldl'
        (\vs (nextSteps,pushTime) ->  applyResults (runSteps nextSteps vs) pushTime vs)
        prev
        -- this third element in this below is a waste every time EXCEPT the last time step.
        $ map (\(x,r) -> (x,r && rerun) ) [(supportCalcSteps,False),(calcSteps,True),(supportCalcSteps,False)] 

runTimeSteps_Print =
    foldM_
        (\prev step -> 
            let next = foldl' (\lastState _ -> runSingleStep lastState False) (runSingleStep prev True) [1..2]
                timeLevel = timePos $ head $ calculatedPositions next -- the first element of the list of calculated positions will give the current time slot
                backOneTimeLevel = pushBackTime timeLevel 
            in do 
                --{-
                -- use the most recent valset but one time step BACK because that will have been corrected. the most recent time level is NOT YET CORRECTED. it has JUST BEEN CALCULATED AND IS INTERMEDIATE
                mapM_
                    (\nextProp ->
                        GP.plotAsync ( PNG.cons $ "c:\\temp\\"++ show nextProp ++ "\\"++show step ++".png") 
                        $! plotDomainLinear $! valSetToGrid next backOneTimeLevel nextProp (1+maxPos Y) (quot (maxPos Z)  2)
                    ) $ enumFrom Speed
                mapM_
                    (\nextProp ->
                        GP.plotAsync ( PNG.cons $ "c:\\temp\\"++ show nextProp ++ " LOGSCALE" ++ "\\"++show step ++".png") 
                        $! plotDomainLog $! valSetToGrid next backOneTimeLevel nextProp (1+maxPos Y) (quot (maxPos Z)  2)
                    ) $ enumFrom Speed
                ---}
                putStrLn $ show $ length (calculatedPositions next)
                putStrLn $ show $ Map.size (vals next)
                putStrLn $ "step: " ++ show step 
                putStrLn $ "timeLevel: " ++ show timeLevel
                putStrLn $ "timeLevelAbsolute: " ++ show (timeLevelAbsolute next)
                putStrLn " "
                --putStrLn $ stringDomain U timeLevel (1 + maxPos Y) prev (quot (maxPos Z)  2)
                return next
        )
        initialGrid
        [0..999999]
 
testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

testEq :: forall a.Reader (ValSet Double) (Equation (Position -> [a])) -> Equation a
testEq eq = getDiscEqInstance ( runReader eq $! initialGrid) testPosition 
            
--writeTermsOrig:: (Num a, Show a, Fractional a)=> [Term a] -> String
writeTermsOrig terms vs=
    let writeTerm prev t= prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) vs ++ " + "
            Derivative d _ side _-> 
                "approxDeriv ( " ++ writeTerms [approximateDerivative t testPosition vs] vs
                    ++ show d ++ " " ++ show  side ++" ) + " 
    in foldl' writeTerm " " (reverse terms )  

writeTerms :: [Term Double] -> ValSet Double-> String
writeTerms terms vs=
    let (_:_:xs) =  reverse $! writeTermsOrig terms vs
    in reverse xs
  
testPosition =   Position [1, 1, 0] 3 0
    
makeRows :: [[Double]] -> [Double]-> [Double] -> Int -> Int-> [[Double]]    
makeRows whole curr [] _ _ = whole ++ [curr]    
makeRows whole curr items 0 width = makeRows (whole ++ [curr] ) [] items width width          
makeRows whole curr (x:xs) r width= makeRows whole (curr++[x]) xs (r-1) width   

valSetToGrid vs timeLevel property rowLength zLevel=
    let positions = 
            filter (\next-> getPositionComponent next Z == zLevel) $ 
            map 
                (\(Position p d _) -> Position p d timeLevel)
                makeAllPositions
    in makeRows [] [] 
            (map (\next -> prop Nondirectional property next Center vs )  positions )
            rowLength
            rowLength 
             
--stringDomain:: (Num a, Fractional a, Show a ) => Property ->[Position]->Int-> ValSet a -> String
stringDomain property timeLevel rowLength set zLevel =
    let rows = valSetToGrid set timeLevel property rowLength zLevel  
        strRows = map ( foldl' (\prev next-> prev ++ " " ++ show next) "" ) rows
    in foldl' (\prev next -> prev ++ "\n" ++ next ) "" strRows 
            
main:: IO()
main = 
    putStrLn "starting ..... "
    >>= (\_-> print ( solveUnknown testEquation ( Position [0, 0, 0] 3 0) initialGrid )) 
    >>= (\_ -> putStrLn $ writeTerms (distributeMultiply testTerms 2) initialGrid)
    >>= (\_ -> print $ prop Nondirectional U testPosition Center  initialGrid)
    >>= (\_ -> putStrLn " continuity ------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq continuity ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq continuity) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq continuity) testPosition initialGrid )
    >>= (\_ -> putStrLn " U Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq uMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq uMomentum ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq uMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " V Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq vMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq vMomentum ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq vMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " W Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq wMomentum ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq wMomentum ) initialGrid ) 
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq wMomentum) testPosition initialGrid )
    >>= (\_ -> putStrLn " ENERGY ------------ ")
    >>= (\_ -> putStrLn $ writeTerms ( rhs $ testEq energy ) initialGrid )
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms ( lhs $ testEq energy ) initialGrid )
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq energy) testPosition initialGrid )
    -- >>= (\_ -> putStrLn $! stringDomain Pressure ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain Pressure (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    >>= (\_ -> runTimeSteps_Print )
