
module Main where

import qualified Data.Map.Strict as Map 
import Data.Maybe
import Control.Monad.Reader as Reader
-- import Control.Monad.Par as Par
import SolutionDomain
import CalculusUtils
import Data.List

-- instance Par.NFData Property 
-- instance Par.NFData Position

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
        ([Constant 2, Constant 3 ] |> addTerms [Unknown (-0.025),Unknown (-0.05)])
        U

continuity::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
continuity = do
    env <- ask
    return $ let integrate = integUnknown env Density 
        in Equation
            ((integrate Temporal [runReader drho_dt env] >>= integrate Spatial) 
              : integrateTerms integrate env ( divergenceWithProps [Density]) )
            [ const [Constant 0]] 
            Density
    
uMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
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

vMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
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

wMomentum::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
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

energy::(Num a, Fractional a, RealFloat a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
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

gasLawPressure::(Num a, Fractional a)=> Reader (ValSet a) (Equation (Position-> [Term a]))
gasLawPressure = do
    env <- ask
    return $  
        Equation
            [ const [Unknown 1]]
            [ \pos -> [ Constant $ gasConstantR * runReader(prop Density pos Center ) env 
                * runReader (prop Temperature pos Center) env ] ]
            Pressure          
    
getDiscEqInstance:: Equation (Position -> [Term a]) -> Position -> Equation (Term a)
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) $! l) (concatMap (\t -> t pos) $! r) up
    
advanceTime :: Position -> Position
advanceTime (Position x y z t ) = Position x y z (t+1)    
    
--applyDiffEq :: (Fractional a, NFData a)=>
applyDiffEq :: (Fractional a)=> 
    ValSet a -> Equation (Position -> [Term a]) -> Bool-> [ (Position,Property,a,Bool)]    
applyDiffEq (ValSet p _ _ _) eq saveAtNextTime=
    -- runPar $ parMap
    map
        (\pos -> 
            let discEquation = getDiscEqInstance eq pos 
                solvedProperty = unknownProperty discEquation
                newValue = solveUnknown discEquation pos
            in ( pos, solvedProperty, newValue, saveAtNextTime)
        ) p
  
applyResults ::  [(Position, Property, a,Bool)]-> ValSet a -> ValSet a
applyResults res (ValSet p v av sl) = 
    let newVals = foldl' 
            (\prev (pos,property,newVal,saveAtNextTime)  ->
                let newPos = if saveAtNextTime 
                        then advanceTime pos 
                        else pos
                    subdict = fromMaybe Map.empty (Map.lookup newPos prev)  
                in Map.insert newPos (Map.insert property newVal subdict) prev
            ) v res
    in ValSet (map advanceTime p) newVals av sl

calcSteps :: (Fractional a, RealFloat a)=> 
    [(Reader (ValSet a) (Equation (Position-> [Term a])) , Bool)]
calcSteps = [ 
    (gasLawPressure, False) 
    ,(continuity, True)
    ,(uMomentum, True)
    ,(vMomentum,  True)
    ,(wMomentum,  True)
    ,(energy,  True)  
    ] 

runSingleStep prev _ = 
    let results = 
            -- runPar $ parMap 
            map
                (\x-> applyDiffEq prev (runReader (fst x) prev) (snd x) ) 
                calcSteps 
    in applyResults (concatMap id results) prev
                
runTimeSteps :: ValSet Double
runTimeSteps = (\x -> foldl'  runSingleStep x [0..0] ) $! initialGrid  
 
testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

testEq:: (Num a, Fractional a)=> Reader (ValSet a ) (Equation (Position ->[Term a])) -> Equation (Term a)
testEq eq = getDiscEqInstance ( runReader eq $! initialGrid) testPosition 
            
writeTermsOrig:: (Num a, Show a, Fractional a)=> [Term a] -> String
writeTermsOrig terms =
    let writeTerm prev t= prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) ++ " + "
            Derivative d _ side _-> 
                "approxDeriv ( " ++ writeTerms [approximateDerivative t testPosition ] 
                    ++ show d ++ " " ++ show  side ++" ) + " 
    in foldl' writeTerm " " (terms |> reverse)  

writeTerms terms =
    let (_:_:xs) = writeTermsOrig terms |> reverse
    in xs |> reverse
  
testPosition =   Position 4 5 0 0
    
makeRows :: [[a]] -> [a]-> [a] -> Int -> Int-> [[a]]    
makeRows whole curr [] _ _ = whole ++ [curr]    
makeRows whole curr items 0 width = makeRows (whole ++ [curr] ) [] items width width          
makeRows whole curr (x:xs) r width= makeRows whole (curr++[x]) xs (r-1) width   
             
stringDomain:: (Num a, Fractional a, Show a ) => Property ->[Position]->Int-> ValSet a -> String
stringDomain property positions rowLength set =
    let rows = 
            makeRows [[]] [] 
                (map (\next -> runReader (prop property next Center) set ) positions )
                rowLength
                rowLength
        strRows = map ( foldl' (\prev next-> prev ++ " " ++ show next) "" ) rows
    in foldl' (\prev next -> prev ++ "\n" ++ next ) "" strRows 
            
main:: IO()
main = 
    putStrLn "starting ..... "
    >>= (\_ -> putStrLn $ stringDomain U (calculatedPositions runTimeSteps) (1+maxPos X) runTimeSteps  )
