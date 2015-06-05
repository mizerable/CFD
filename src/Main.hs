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

timeStep::Double
timeStep = 0.0001 

maxPos:: Direction -> Int
maxPos d = case d of 
    X -> 100000
    Y -> 100000
    Z -> 100000

getPositionComponent (Position x y z _) d = case d of 
    X -> x
    Y -> y
    Z -> z

modifyPositionComponent (Position x y z t) direction amt= case direction of 
    X -> Position (x+amt) y z t
    Y -> Position x (y+amt) z t 
    Z -> Position x y (z+amt) t
    
offsetPosition:: Position->Side ->Position
offsetPosition (Position x y z t) side = case side of
    Center -> Position x y z t 
    Now -> Position x y z t 
    Prev -> Position x y z (max 0 (t - timeStep))  
    _ -> 
        let maxOrMin = if isUpperSide side then min else max
            offsetAmount = if isUpperSide side then 1 else (-1)
            direction = directionFromCenter side
            boundary = if isUpperSide side 
                then maxPos direction 
                else 0
        in modifyPositionComponent (Position x y z t) direction offsetAmount   
        
    --East -> Position (min (maxPos X) (x+1)) y z t
    --West -> Position (max 0 (x-1)) y z t
    --Top -> Position x y (min(maxPos Z) (z+1)) t
    --Bottom -> Position x y (max 0 (z-1)) t
    --North -> Position x (min (maxPos Y) (y+1)) z t
    --South -> Position x (max 0 (y-1)) z t

prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
prop property position side = do
    (ValSet _ v _ _) <- ask 
    return $ 
        let neighbor = offsetPosition position side
            getVal p = fromJust $ Map.lookup p v >>= Map.lookup property  
        in average [getVal position,getVal neighbor]

initialGrid::(Num a) => ValSet a
initialGrid= undefined
         
direcDimenType:: Direction -> DimensionType
direcDimenType direc = case direc of
    Time -> Temporal
    _ -> Spatial

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = length terms |> fromIntegral
        f n p = p + n / len
    in foldr f 0 terms

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

volumeOrInterval:: DimensionType -> Position -> Reader (ValSet Double) Double
volumeOrInterval dimetype position = do 
    vs <- ask
    return $ case dimetype of
        Temporal -> timeStep
        Spatial -> enumFrom X |> foldr (\d -> \p -> p * runReader (sideLength d position) vs ) 1

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
            sl =  runReader (sideLength direction position)
            sln = runReader (sideLength direction neighbor)
            interval = average [sl vs, sln vs]
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
                    (1/ sl vs)    
    _ -> []

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

integSurface:: (Num a)=> (Side->a) -> Position -> Direction -> Reader (ValSet a) [Term a]
integSurface f position direction = do
    vs <- ask
    return $ 
        let sides = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f else (\x-> x * (-1)).f 
                in (case (direcDimenType direction,isUpper) of
                    (Temporal,True) -> Unknown
                    _ -> Constant)
                    ( modf s*( runReader (sideArea s position) vs))
        in [value (fst sides) True , value (snd sides) False]       
       
integSingleTerm::  Term Double -> DimensionType -> Position -> Reader (ValSet Double) [Term Double]
integSingleTerm term dimetype cellposition =  do
    vs <- ask
    return $ case term of
        Derivative direction func _ ->  
            if dimetype == direcDimenType direction 
                then runReader (integSurface ( func cellposition ) cellposition direction) vs
                else distributeMultiply [term] $ runReader (volumeOrInterval dimetype cellposition) vs 
        _ -> distributeMultiply [term] $ runReader (volumeOrInterval dimetype cellposition) vs  

integ::  ValSet Double -> DimensionType -> [Term Double] ->Position -> [Term Double]
integ vs dimetype terms cellposition = case terms of
    [] -> []
    (x:xs) -> runReader (integSingleTerm x dimetype cellposition) vs ++ integ vs dimetype xs cellposition   
      
drho_dt = do
    d<-ask 
    return $ Derivative Time (\x -> \s-> runReader (prop Density x s) d ) Center

drhodu_dt =do
    d<-ask 
    return $ Derivative Time (\x-> \s -> runReader (prop Density x s) d * runReader (prop U x s) d) Center
 
drhodv_dt = do
    d<-ask   
    return $ Derivative Time (\x-> \s -> runReader (prop Density x s) d* runReader (prop V x s) d) Center

drhodw_dt = do
    d<-ask  
    return $ Derivative Time (\x-> \s -> runReader (prop Density x s) d* runReader (prop W x s) d) Center  

-- the only thing returned is the new value for the density at this position 
continuity:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            [ integrate Temporal [runReader drho_dt env] >>= integrate Spatial ,
              integrate Temporal [runReader drho_dt env] >>= integrate Spatial, 
              integrate Spatial [runReader drhodu_dt env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodv_dt env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodw_dt env] >>= integrate Temporal ] 
            [\_ -> [Constant 0]] 
    
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
    

