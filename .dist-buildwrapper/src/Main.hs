module Main where

import qualified Data.Map as Map 
import Data.Maybe

data Side = East | West | North | South | Top | Bottom | Now | Prev | Center deriving (Show,Eq,Ord)
data Direction = Time | X | Y | Z deriving (Enum)
data DimensionType = Temporal | Spatial
data Equation a = Equation{
    rhs::[Term a]
    ,lhs::[Term a] }
data IntegralType = Body | Surface deriving (Eq)
data Property = U | V | W | Density | Temperature
data Position = Position {
    xPos::Int
    ,yPos::Int
    ,zPos::Int
    ,timePos::Int } deriving (Eq, Ord)
data Derivative a = Derivative{
    denom::Direction
    ,numer::Property
    ,function:: Position->Side->a }
data ValSet a = ValSet{
    vals:: Map.Map Position (Map.Map Property a)
    ,areaVal::Map.Map Position (Map.Map Side a) }
data Expression a = Expression{getTerms::[Term a]}
data Term a = Constant {val::a} | Unknown { coeff::a } | SubExpression {expression::Expression a} 

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

solveUnknown::(Fractional a)=> Equation a->a
solveUnknown equation = 
    let sumUnknown n p =  p + case n of
            Unknown _-> coeff n
            SubExpression _ -> sumExpression sumUnknown (expression n |> getTerms)
            _ -> 0
        sumConstants n p =  p + case n of
            Constant _-> val n
            SubExpression _ -> sumExpression sumConstants (expression n |> getTerms)
            _ -> 0
        sumExpression s e = foldr s 0 e
        lhsUnknown = sumExpression sumUnknown (lhs equation)
        rhsUnknown = sumExpression sumUnknown (rhs equation)
        lhsConstants = sumExpression sumConstants (lhs equation)
        rhsConstants = sumExpression sumConstants (rhs equation)
    in (rhsConstants - lhsConstants)/(lhsUnknown-rhsUnknown)
        
testEquation:: Equation Double
testEquation = 
    Equation 
        --[Unknown 2, Constant 4, Constant (-0.232)]
        [Unknown 2, SubExpression (Expression [Constant 2, Constant 2, Unknown 0.025]), Constant (-0.232)] 
        ([Constant 2, Constant 3 ] |> addTerms [Unknown (-0.025),Unknown (-0.05)])

direcIntegType:: Direction -> IntegralType
direcIntegType direction = case direction of 
    Time -> Body
    _ -> Surface

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
            SubExpression _ -> distributeMultiply ( (getTerms.expression) term  ) m
    in undefined    

prop::Property->Side->Position->a
prop property side position = undefined

integSurface:: (Num a)=> (Side->a) -> Position -> Direction -> [Term a]
integSurface f position direction =
    let sides = boundaryPair direction 
        value s = Constant (f s*(sideArea s position |> fromJust))
    in [value (fst sides) , value (snd sides)]       
       
integ:: Derivative Double -> Direction -> Position ->[Term Double]
integ derivative direction cellposition
    | direcIntegType direction == direcIntegType (denom derivative)
        = integSurface ( cellposition |> function derivative) cellposition direction 
    | otherwise =  [Constant (function derivative cellposition Center * 
        volumeOrInterval (direcDimenType direction) cellposition )]            
       
-- the only thing returned is the new value for the density at this position 
continuity::(Num a)=> Position -> a
continuity pos =
    -- drho/dt 
    undefined

applyContinuity::(Num a)=>ValSet a->ValSet a 
applyContinuity valset = undefined


main:: IO()
main = putStrLn ( (show.solveUnknown) testEquation )

