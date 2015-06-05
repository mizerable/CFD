module Main where

import qualified Data.Map as Map 
import Data.Maybe
import Control.Monad.State as State
import Control.Monad.Reader as Reader
  
data Side = East | West | North | South | Top | Bottom | Now | Prev | Center deriving (Show,Eq,Ord)
data Direction = Time | X | Y | Z deriving (Enum,Ord,Eq)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    rhs::[a]
    ,lhs::[a]}    
data IntegralType = Body | Surface deriving (Eq)
data Property = U | V | W | Density | Temperature deriving (Ord,Eq)
data Position = Position {
    xPos::Int
    ,yPos::Int
    ,zPos::Int
    ,timePos::Double } deriving (Eq, Ord)
data ValSet a = ValSet{
    positions:: [Position]
    ,vals:: Map.Map Position (Map.Map Property a)
    ,areaVal::Map.Map Position (Map.Map Side a)
    ,sideLen:: Map.Map Position (Map.Map Direction a) }
data Expression a = Expression{getTerms::[Term a]} 
data Term a = Constant {val::a} | Unknown { coeff::a } | SubExpression {expression::Expression a} 
    | Derivative { denom::Direction ,function:: Position->Side->a, centered::Side }

direcDimenType:: Direction -> DimensionType
direcDimenType direc = case direc of
    Time -> Temporal
    _ -> Spatial

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = length terms |> fromIntegral
        f n p = p + n / len
    in foldr f 0 terms

timeStep::Double
timeStep = 0.0001 

isUpperSide:: Side -> Bool
isUpperSide side = case side of
    East -> True
    North -> True
    Top -> True
    Now -> True 
    _->False

boundaryPair:: Direction -> (Side,Side)
boundaryPair d = case d of 
     X -> (East,West)
     Y -> (North,South)
     Z -> (Top,Bottom)
     Time -> (Now,Prev)

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

directionFromCenter:: Side-> Direction
directionFromCenter side = case side of
    East -> X
    West -> X
    North -> Y
    South -> Y
    Top -> Z
    Bottom -> Z
    Now -> Time
    Prev -> Time

volumeOrInterval:: ValSet Double -> DimensionType -> Position -> Double
volumeOrInterval vs dimetype position = case dimetype of
    Temporal -> timeStep
    Spatial -> enumFrom X |> foldr (\d -> \p -> p * sideLength vs d position |> fromJust ) 1

addTerm:: [Term a]->Term a->[Term a]
addTerm terms term = terms ++ [term]

addTerms:: [Term a]-> [Term a]-> [Term a]
addTerms terms1 terms2 = case terms2 of
    ([x]) -> addTerm terms1 x
    (x:xs) -> addTerms (addTerm terms1 x) xs
    _ -> terms1

approximateDerivative::(Num a, Fractional a)=> ValSet a -> Term a -> Position-> [Term a]
approximateDerivative vs deriv position= case deriv of 
    (Derivative direction func side) ->
        let neighbor = offsetPosition position side 
            interval = average [ sideLength vs direction position |> fromJust, sideLength vs direction neighbor |> fromJust]
            thisVal = func position Center
            neighborVal = func neighbor Center
            f first = if isUpperSide side && first then neighborVal else thisVal
        in if direction == directionFromCenter side 
            then distributeMultiply [ Constant (f True) , Constant ( (-1)*f False) ] (1/interval)
            else 
                let proxyNeighbor p t = func (offsetPosition p $ t $ boundaryPair direction) Center 
                    thisProxyNeighborA = proxyNeighbor position fst
                    thisProxyNeighborB = proxyNeighbor position snd
                    neighborProxyNeighborA = proxyNeighbor neighbor fst
                    neighborProxyNeighborB = proxyNeighbor neighbor snd 
                in distributeMultiply [ Constant (average [thisProxyNeighborA,neighborProxyNeighborA])
                    , Constant ((-1)* average [thisProxyNeighborB,neighborProxyNeighborB]) ] 
                    (1/sideLength vs direction position|> fromJust)    
    _ -> []

offsetPosition:: Position->Side ->Position
offsetPosition position side = undefined

solveUnknown::(Fractional a)=> ValSet a->Equation (Term a)->Position->a
solveUnknown vs (Equation l r) position= 
    let sumUnknown n p =  p + case n of
            Unknown u-> u
            SubExpression s -> sumExpression sumUnknown $ getTerms s
            _ -> 0
        sumConstants n p =  p + case n of
            Constant c-> c
            Derivative _ _ _-> sumExpression sumConstants $ approximateDerivative vs n position  
            SubExpression s -> sumExpression sumConstants $ getTerms s
            _ -> 0
        sumExpression s e = foldr s 0 e
        lhsUnknown = sumExpression sumUnknown l
        rhsUnknown = sumExpression sumUnknown r
        lhsConstants = sumExpression sumConstants l
        rhsConstants = sumExpression sumConstants r
    in (rhsConstants - lhsConstants)/(lhsUnknown-rhsUnknown)
        
testEquation:: Equation (Term Double)
testEquation = 
    Equation 
        --[Unknown 2, Constant 4, Constant (-0.232)]
        [Unknown 2, SubExpression (Expression [Constant 2, Constant 2, Unknown 0.025]), Constant (-0.232)] 
        ([Constant 2, Constant 3 ] |> addTerms [Unknown (-0.025),Unknown (-0.05)])

(|>)::a->(a->b)->b
(|>) x y = y x

sideArea:: Side -> Position -> Reader (ValSet a) a
sideArea s position = do
    vs <- ask  
    return $ fromJust $ areaVal vs |> Map.lookup position >>= Map.lookup s   

sideLength:: Direction -> Position -> Reader (ValSet a) a
sideLength d position = do
    vs <- ask
    return $ fromJust $ sideLen vs |> Map.lookup position >>= Map.lookup d

prop:: Property->Position->Side-> Reader (ValSet a) a
prop property position side = undefined   

initialGrid::(Num a) => ValSet a
initialGrid= undefined

distributeMultiply::(Num a)=> [Term a]->a->[Term a]
distributeMultiply terms m =
    let mult term = case term of
            Constant c -> [Constant (c * m)]
            Unknown u -> [Unknown (u * m)]
            SubExpression s -> distributeMultiply ( getTerms s  ) m
            Derivative direc func side-> 
                let modf x s =  func x s * m
                in [Derivative direc modf side]
    in concatMap mult terms    

integSurface:: (Num a)=> ValSet a -> (Side->a) -> Position -> Direction -> [Term a]
integSurface vs f position direction =
    let sides = boundaryPair direction 
        value s isUpper =
            let modf = if isUpper then f else (\x-> x * (-1)).f 
            in (case (direcDimenType direction,isUpper) of
                (Temporal,True) -> Unknown
                _ -> Constant)
                ( modf s*(sideArea vs s position |> fromJust))
    in [value (fst sides) True , value (snd sides) False]       
       
integSingleTerm:: ValSet Double -> Term Double -> DimensionType -> Position ->[Term Double]
integSingleTerm vs term dimetype cellposition =  case term of
    Derivative direction func _ ->  
        if dimetype == direcDimenType direction then integSurface vs ( func cellposition ) cellposition direction
                else distributeMultiply [term] $ volumeOrInterval vs dimetype cellposition   
    _ -> distributeMultiply [term] $ volumeOrInterval vs dimetype cellposition 

integ::  ValSet Double -> DimensionType -> [Term Double] ->Position -> [Term Double]
integ vs dimetype terms cellposition = case terms of
    [] -> []
    (x:xs) -> integSingleTerm vs x dimetype cellposition ++ integ vs dimetype xs cellposition   
       
drho_dt:: (Num a)=> ValSet a-> Term a       
drho_dt d =  Derivative Time (prop d Density) Center

drhodu_dt:: (Num a)=> ValSet a->Term a
drhodu_dt d =  Derivative Time (\x-> \s -> prop d Density x s* prop d U x s) Center

drhodv_dt:: (Num a)=> ValSet a->Term a 
drhodv_dt d =  Derivative Time (\x-> \s -> prop d Density x s*prop d V x s) Center

drhodw_dt:: (Num a)=> ValSet a->Term a
drhodw_dt d =  Derivative Time (\x-> \s -> prop d Density x s*prop d W x s) Center  

-- the only thing returned is the new value for the density at this position 
continuity:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = asks
    (\env -> 
        let integrate = integ env 
        in Equation
            [ integrate Temporal [drho_dt env] >>= integrate Spatial ,
              integrate Temporal [drho_dt env] >>= integrate Spatial, 
              integrate Spatial [drhodu_dt env] >>= integrate Temporal, 
              integrate Spatial [drhodv_dt env] >>= integrate Temporal, 
              integrate Spatial [drhodw_dt env] >>= integrate Temporal ] 
            [\_ -> [Constant 0]] )
    
uMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
uMomentum = undefined

vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
vMomentum = undefined    

wMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
wMomentum = undefined    

energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
energy = undefined    

gasLaw:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
gasLaw= undefined        
    
getUnknownPropertyType:: Equation a -> Property
getUnknownPropertyType = undefined    
    
getDiscEqInstance:: Equation (Position -> [Term a]) -> Position -> Equation (Term a)
getDiscEqInstance (Equation l r) pos = Equation (concatMap (\t -> t pos) l) (concatMap (\t -> t pos) r)
    
applyDiffEq :: (Fractional a)=>ValSet a -> Equation (Position -> [Term a]) -> ValSet a    
applyDiffEq (ValSet p v av sl) eq =
    let newVals = foldr
            (\pos -> \dict -> 
                let subDict = fromJust $ Map.lookup pos dict
                    discEquation=getDiscEqInstance eq pos 
                    solvedProperty = getUnknownPropertyType discEquation
                    newValue = solveUnknown (ValSet p v av sl) discEquation pos  
                in Map.insert pos (Map.insert solvedProperty newValue subDict)dict )  
            v p 
    in ValSet p newVals av sl
    
updateDomain::(Fractional a)=> Reader (ValSet a) (Equation (Position -> [Term a])) -> State (ValSet a) ()
updateDomain equation = state $ \prev -> ((),applyDiffEq prev $ runReader equation prev)       
  
runTimeSteps:: State (ValSet Double) [()]
runTimeSteps =  mapM 
        (\_ ->  updateDomain continuity  
            >>= \_ -> updateDomain uMomentum
            >>= \_ -> updateDomain vMomentum
            >>= \_ -> updateDomain wMomentum
            >>= \_ -> updateDomain energy
            >>= \_ -> updateDomain gasLaw ) 
        [1..10] 
    
testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

writeTerms:: (Num a, Show a)=> [Term a] -> String
writeTerms terms =
    let writeTerm t prev = prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) ++ " + "
    in foldr writeTerm " " (terms |> reverse)  

main:: IO()
main = 
    print ( solveUnknown initialGrid testEquation $ Position 0 0 0 0.0) 
    >>= \_ -> putStrLn $ writeTerms $ distributeMultiply testTerms 2
    

