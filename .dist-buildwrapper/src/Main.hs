module Main where

import qualified Data.Map as Map 
import Data.Maybe

data Side = East | West | North | South | Top | Bottom | Now | Prev | Center deriving (Show,Eq,Ord)
data Direction = Time | X | Y | Z deriving (Enum)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    rhs::[a]
    ,lhs::[a]}    
data IntegralType = Body | Surface deriving (Eq)
data Property = U | V | W | Density | Temperature
data Position = Position {
    xPos::Int
    ,yPos::Int
    ,zPos::Int
    ,timePos::Double } deriving (Eq, Ord)
data ValSet a = ValSet{
    vals:: Map.Map Position (Map.Map Property a)
    ,areaVal::Map.Map Position (Map.Map Side a) }
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

boundaryPair:: Direction -> (Side,Side)
boundaryPair d = case d of 
     X -> (East,West)
     Y -> (North,South)
     Z -> (Top,Bottom)
     Time -> (Now,Prev)

volumeOrInterval:: DimensionType -> Position -> Double
volumeOrInterval dimetype position = case dimetype of
    Temporal -> timeStep
    Spatial ->
        let getSide f x = sideArea ((f.boundaryPair) x) position |> fromJust 
        in enumFrom X
            |> map (\x-> average [ getSide fst x , getSide snd x])
            |> average

addTerm:: [Term a]->Term a->[Term a]
addTerm terms term = terms ++ [term]

addTerms:: [Term a]-> [Term a]-> [Term a]
addTerms terms1 terms2 = case terms2 of
    ([x]) -> addTerm terms1 x
    (x:xs) -> addTerms (addTerm terms1 x) xs
    _ -> terms1

getSubExpression = getTerms.expression

approximateDerivative:: Term a -> [Term a]
approximateDerivative = undefined

solveUnknown::(Fractional a)=> Equation (Term a)->a
solveUnknown equation = 
    let sumUnknown n p =  p + case n of
            Unknown _-> coeff n
            SubExpression _ -> sumExpression sumUnknown $ getSubExpression n
            _ -> 0
        sumConstants n p =  p + case n of
            Constant _-> val n
            Derivative _ _ _-> sumExpression sumConstants $ approximateDerivative n  
            SubExpression _ -> sumExpression sumConstants $ getSubExpression n
            _ -> 0
        sumExpression s e = foldr s 0 e
        lhsUnknown = sumExpression sumUnknown (lhs equation)
        rhsUnknown = sumExpression sumUnknown (rhs equation)
        lhsConstants = sumExpression sumConstants (lhs equation)
        rhsConstants = sumExpression sumConstants (rhs equation)
    in (rhsConstants - lhsConstants)/(lhsUnknown-rhsUnknown)
        
testEquation:: Equation (Term Double)
testEquation = 
    Equation 
        --[Unknown 2, Constant 4, Constant (-0.232)]
        [Unknown 2, SubExpression (Expression [Constant 2, Constant 2, Unknown 0.025]), Constant (-0.232)] 
        ([Constant 2, Constant 3 ] |> addTerms [Unknown (-0.025),Unknown (-0.05)])

(|>)::a->(a->b)->b
(|>) x y = y x

sideArea s position =  areaVal grid |> Map.lookup position >>= Map.lookup s

grid::(Num a) => ValSet a
grid = undefined

distributeMultiply::(Num a)=> [Term a]->a->[Term a]
distributeMultiply terms m =
    let mult term = case term of
            Constant _ -> [Constant (val term * m)]
            Unknown _ -> [Unknown (coeff term * m)]
            SubExpression _ -> distributeMultiply ( getSubExpression term  ) m
            Derivative _ _ _-> 
                let modf x s =  function term x s * m
                in [Derivative (denom term) modf (centered term)]
    in concatMap mult terms    

prop::Property->Position->Side->a
prop property side position = undefined

integSurface:: (Num a)=> (Side->a) -> Position -> Direction -> [Term a]
integSurface f position direction =
    let sides = boundaryPair direction 
        value s isUpper =
            let modf = if isUpper then f else (\x-> x * (-1)).f 
            in (case (direcDimenType direction,isUpper) of
                (Temporal,True) -> Unknown
                _ -> Constant)
                ( modf s*(sideArea s position |> fromJust))
    in [value (fst sides) True , value (snd sides) False]       
       
integSingleTerm:: Term Double -> DimensionType -> Position ->[Term Double]
integSingleTerm term dimetype cellposition =  case term of
    Derivative _ _ _ ->  
        let direction = denom term
        in if dimetype == direcDimenType direction then integSurface ( cellposition |> function term) cellposition direction
            else distributeMultiply [term] $ volumeOrInterval dimetype cellposition   
    _ -> distributeMultiply [term] $ volumeOrInterval dimetype cellposition 

integ::  DimensionType -> [Term Double] ->Position -> [Term Double]
integ dimetype terms cellposition = case terms of
    [] -> []
    (x:xs) -> integSingleTerm x dimetype cellposition ++ integ dimetype xs cellposition   
       
drho_dt:: (Num a)=> Term a       
drho_dt =  Derivative Time (prop Density) Center

drhodu_dt:: (Num a)=> Term a
drhodu_dt = Derivative Time (\x-> \s -> prop Density x s* prop U x s) Center

drhodv_dt:: (Num a)=> Term a 
drhodv_dt = Derivative Time (\x-> \s -> prop Density x s*prop V x s) Center

drhodw_dt:: (Num a)=> Term a
drhodw_dt = Derivative Time (\x-> \s -> prop Density x s*prop W x s) Center  

(>*>):: (c->a)->(a->c->a)->(c->a)
(>*>) prev next = \input -> next (prev input) input  
     
-- the only thing returned is the new value for the density at this position 
continuity:: Equation (Position-> [Term Double])
continuity = Equation
    [ integ Temporal [drho_dt]  >*> (integ Spatial), 
     integ  Spatial [drhodu_dt] >*> (integ Temporal), 
     integ Spatial [drhodv_dt] >*> (integ Temporal), 
     integ  Spatial [drhodw_dt]>*> (integ Temporal) ] 
    []
    
applyContinuity::(Num a)=>ValSet a->ValSet a 
applyContinuity valset = undefined

testTerms = [Unknown 2.4, Constant 1.2, Constant 3.112, Unknown (-0.21),  SubExpression (Expression [Constant 2, Constant 2, SubExpression (Expression [Unknown 0.33333])])]

writeTerms:: (Num a, Show a)=> [Term a] -> String
writeTerms terms =
    let writeTerm t prev = prev ++ case t of
            Unknown _ -> show (coeff t) ++ "X + "
            Constant _ -> show (val t) ++ " + "
            SubExpression _ -> writeTerms (getSubExpression t) ++ " + "
    in foldr writeTerm " " (terms |> reverse)  

main:: IO()
main = 
    putStrLn ( (show.solveUnknown) testEquation )
    >>= \_ -> putStrLn $ writeTerms $ distributeMultiply testTerms 2


