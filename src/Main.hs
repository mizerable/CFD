
module Main where

import qualified Data.Map.Strict as Map 
import Data.Maybe
import Control.Monad.Reader as Reader
import Control.Monad.Par as Par
import SolutionDomain
import CalculusUtils
import GeometryStuff
import Data.List

import qualified Graphics.Gnuplot.Advanced as GP
import qualified Graphics.Gnuplot.Plot.TwoDimensional as Plot2D
import qualified Graphics.Gnuplot.Graph.TwoDimensional as Graph2D
import Graphics.Gnuplot.Plot.TwoDimensional (linearScale, )
import Graphics.Gnuplot.Terminal.PNG as PNG 

instance Par.NFData Property 
instance Par.NFData Position

defaultInflow:: (Num a, Fractional a) => a
defaultInflow = 1

gasConstantR :: (Num a, Fractional a) => a
gasConstantR = 8.314

specificHeatCp :: (Num a, Fractional a ) => a
specificHeatCp = gasConstantR + specificHeatCv

heatConductivityK:: (Num a, Fractional a ) => a
heatConductivityK = 0.1

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
continuity::Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = do
    env <- ask
    return $ let integrate = integUnknown env Density 
        in Equation
            ((integrate Temporal [runReader drho_dt env] >>= integrate Spatial) 
              : integrateTerms integrate env ( divergenceWithProps [Density]) )
            [ const [Constant 0]] 
            Density
    
--uMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
uMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
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

vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
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

wMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
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

energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
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

gasLawPressure::Reader (ValSet Double) (Equation (Position-> [Term Double]))
gasLawPressure = do
    env <- ask
    return $  
        Equation
            [ const [Unknown 1]]
            [ \pos -> [ Constant $ gasConstantR * prop Density pos Center  env 
                * prop Temperature pos Center env ] ]
            Pressure          
    
getDiscEqInstance:: Equation (Position -> [Term a]) -> Position -> Equation (Term a)
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) $! l) (concatMap (\t -> t pos) $! r) up
    
advanceTime :: Position -> Position
advanceTime (Position x y z t ) = Position x y z $! mod (t+1) storedSteps     
    
--applyDiffEq :: (Fractional a, NFData a)=> ValSet a -> Equation (Position -> [Term a]) -> Bool-> [ (Position,Property,a,Bool)]    
applyDiffEq (eq, saveAtNextTime,getPos) env =
    runPar $ parMap
    --map
        (\pos -> 
            let discEquation = getDiscEqInstance (runReader eq env) pos 
                solvedProperty = unknownProperty discEquation
                newValue = solveUnknown discEquation pos
            in ( pos, solvedProperty, newValue, saveAtNextTime)
        ) (getPos env)
  
-- applyResults ::  [(Position, Property, a,Bool)]-> ValSet a -> ValSet a
applyResults res pushTime (ValSet p v av sl) = 
    let newVals = foldl' 
            (\prev (pos,property,newVal,saveAtNextTime)  ->
                let newPos = if saveAtNextTime 
                        then advanceTime pos 
                        else pos
                    subdict = fromMaybe Map.empty (Map.lookup newPos prev)  
                in Map.insert newPos (Map.insert property newVal subdict) prev
            ) v res
    in ValSet 
        (if pushTime then map advanceTime p else p) 
        newVals av sl

-- calcSteps :: (Fractional a, RealFloat a)=>  [(Reader (ValSet a) (Equation (Position-> [Term a])) , Bool)]
calcSteps = [ 
    (continuity, True, allPositionsCurrentTime )
    ,(uMomentum, True, calculatedPositions )
    ,(vMomentum,  True, calculatedPositions )
    ,(wMomentum,  True, calculatedPositions )
    ,(energy,  True, allPositionsCurrentTime )  
    ] 

supportCalcSteps = [
    (gasLawPressure, False , allPositionsCurrentTime  )
    ]

allPositionsCurrentTime env = 
    let curTimeLevel = timePos $ head $ calculatedPositions env
    in map (\x-> modifyPositionComponent x Time curTimeLevel ) makeAllPositions 

runSingleStep prev _ = 
    let runSteps steps vs= concat $!
            runPar $ parMap
            --map 
                (\step-> applyDiffEq step vs ) 
                steps  
    in foldl'
        (\vs (nextSteps,pushTime) ->  applyResults (runSteps nextSteps vs) pushTime vs)
        prev
        -- this third element in this below is a waste every time EXCEPT the last time step.
        [(supportCalcSteps,False),(calcSteps,True),(supportCalcSteps,False)] 
                
runTimeSteps :: ValSet Double
runTimeSteps = (\x -> foldl'  runSingleStep x [0..25] ) $! initialGrid  
 
plotDomain :: [[Double]]-> Plot2D.T Int Int
plotDomain grid = 
    Plot2D.function Graph2D.image 
        (liftM2 (,) [0..length grid -1 ] [0.. length (head grid) -1 ]) 
            $ \(x,y) -> (grid!!x)!!y
             
runTimeSteps_Print =
    foldM_
        (\prev step -> do
            let timeLevel = timePos $ head $ calculatedPositions prev
            GP.plotAsync ( PNG.cons $ "c:\\temp\\"++ (show step) ++".png") 
                $! plotDomain $! valSetToGrid prev timeLevel U (1+maxPos X)
            putStrLn $ show $ length (calculatedPositions prev)
            putStrLn $ "step: " ++ show step 
            putStrLn $ "timeLevel: " ++ show timeLevel
            putStrLn $ stringDomain U timeLevel (1 + maxPos X) prev
            return $! runSingleStep prev ()
        )
        initialGrid
        [0..50]
 
testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

--testEq:: (Num a, Fractional a)=> Reader (ValSet a ) (Equation (Position ->[Term a])) -> Equation (Term a)
testEq eq = getDiscEqInstance ( runReader eq $! initialGrid) testPosition 
            
--writeTermsOrig:: (Num a, Show a, Fractional a)=> [Term a] -> String
writeTermsOrig terms =
    let writeTerm prev t= prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) ++ " + "
            Derivative d _ side _-> 
                "approxDeriv ( " ++ writeTerms [approximateDerivative t testPosition ] 
                    ++ show d ++ " " ++ show  side ++" ) + " 
    in foldl' writeTerm " " (reverse terms )  

writeTerms terms =
    let (_:_:xs) =  reverse $! writeTermsOrig terms 
    in reverse xs
  
testPosition =   Position 1 1 0 0
    
makeRows :: [[a]] -> [a]-> [a] -> Int -> Int-> [[a]]    
makeRows whole curr [] _ _ = whole ++ [curr]    
makeRows whole curr items 0 width = makeRows (whole ++ [curr] ) [] items width width          
makeRows whole curr (x:xs) r width= makeRows whole (curr++[x]) xs (r-1) width   

valSetToGrid vs timeLevel property rowLength =
    let positions = map 
            (\(Position x y z _) -> Position x y z timeLevel)
            makeAllPositions
    in makeRows [] [] 
            (map (\next -> prop property next Center vs )  positions )
            rowLength
            rowLength 
             
--stringDomain:: (Num a, Fractional a, Show a ) => Property ->[Position]->Int-> ValSet a -> String
stringDomain property timeLevel rowLength set =
    let rows = valSetToGrid set timeLevel property rowLength  
        strRows = map ( foldl' (\prev next-> prev ++ " " ++ show next) "" ) rows
    in foldl' (\prev next -> prev ++ "\n" ++ next ) "" strRows 
            
main:: IO()
main = 
    putStrLn "starting ..... "
    >>= (\_-> print ( solveUnknown testEquation $ Position 0 0 0 0)) 
    >>= (\_ -> putStrLn $ writeTerms $ distributeMultiply testTerms 2)
    >>= (\_ -> print $ prop U testPosition Center  initialGrid)
    >>= (\_ -> putStrLn " continuity ------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq continuity)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq continuity)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq continuity) testPosition)
    >>= (\_ -> putStrLn " U Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq uMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq uMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq uMomentum) testPosition)
    >>= (\_ -> putStrLn " V Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq vMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq vMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq vMomentum) testPosition)
    >>= (\_ -> putStrLn " W Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq wMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq wMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq wMomentum) testPosition)
    >>= (\_ -> putStrLn " ENERGY ------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq energy)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq energy)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq energy) testPosition)
    >>= (\_ -> putStrLn " Pressure ------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq gasLawPressure)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq gasLawPressure)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> print $ solveUnknown (testEq gasLawPressure) testPosition)
    -- >>= (\_ -> putStrLn $! stringDomain Pressure ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain Pressure (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U ( timePos $ offsetPosition (head $ calculatedPositions runTimeSteps) Prev) (1+maxPos X) runTimeSteps  )
    -- >>= (\_ -> putStrLn $! stringDomain U (timePos $ head $ calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
    >>= (\_ -> runTimeSteps_Print )
