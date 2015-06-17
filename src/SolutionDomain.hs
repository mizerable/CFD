
module SolutionDomain where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Data.Maybe
import Control.Monad.Reader as Reader
import Data.List

data Side =  Now | Prev |East | West | North | South | Top | Bottom | Center deriving (Show,Eq,Ord, Enum)
  
data Direction = Time | X | Y | Z deriving (Enum,Ord,Eq,Show)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    lhs:: ![a]
    ,rhs:: ![a]
    ,unknownProperty:: !Property}    
data IntegralType = Body | Surface deriving (Eq)
data Property = U | V | W | Density | Temperature | Mew | Pressure deriving (Ord,Eq,Enum,Show)
data Position = Position {
    xPos:: !Int
    ,yPos:: !Int
    ,zPos:: !Int
    ,timePos:: !Int } deriving (Eq, Ord, Show)
data ValSet a = ValSet{
    calculatedPositions:: ![Position]
    ,vals:: !(Map.Map Position (Map.Map Property a))
    ,areaVal:: !(Map.Map Position (Map.Map Side a))
    ,sideLen:: !(Map.Map Position (Map.Map Direction a)) }
data Expression a = Expression{getTerms:: ![Term a]} 
data Term a = 
    Constant {val:: !a}
    | Unknown { coeff:: !a }
    | SubExpression {expression:: !(Expression a) }  
    | Derivative { denom:: !Direction ,function:: !(Position->Side->Term a), centered:: !Side, 
        multiplier:: !(Position->Side-> a) } 
        
instance Functor Term  where
    fmap f x = case x of
         Constant c -> Constant $ f c
         Unknown u -> Unknown $ f u
         _ -> undefined    

storedSteps:: Int
storedSteps = 4

maxPos:: Direction -> Int
maxPos d = case d of 
    X -> 15
    Y -> 15
    Z -> 0
    Time -> undefined

removeItems  :: (Ord a, Eq a)=> [a] -> [a]-> [a]
removeItems orig remove= 
    let removeSet = Set.fromList remove
    in filter ( `Set.notMember` removeSet) orig      
           
wallPositionsVals :: ValSet Double
wallPositionsVals = ValSet 
     (obstacle : inflowPositions)
     Map.empty 
     Map.empty Map.empty 

wallPositions :: [Position]
wallPositions = calculatedPositions wallPositionsVals

wallPositionsSet :: Set.Set Position
wallPositionsSet = Set.fromList $! wallPositions 

obstacle = Position (quot (maxPos X)  3) (quot (maxPos Y) 2) 0 0

obstacles = [obstacle]

initialGrid:: ValSet Double
initialGrid= 
    let p = makeAllPositions
        vMap = foldl' (\prev next -> Map.insert next 
            (case next of 
                U-> 1
                V-> 0
                W-> 0
                Density -> 1
                _-> 1
                ) 
            prev) Map.empty (enumFrom U)
        avMap = foldl' (\prev next -> Map.insert next 1 $! prev) Map.empty $!  enumFrom East
        slMap = foldl' (\prev next -> Map.insert next 1 $! prev) Map.empty $! enumFrom X
        v = foldl' (\prev next -> Map.insert next vMap $! prev) Map.empty  $!  p
        av = foldl' (\prev next -> Map.insert next avMap $! prev) Map.empty $!  p
        sl = foldl' (\prev next -> Map.insert next slMap $! prev) Map.empty $!  p
        calcPos = removeItems p $! wallPositions
    in foldl'
        (\prev next -> setVal prev next U 0.0)
        (ValSet calcPos v av sl)
        obstacles
        
setVal:: ValSet a -> Position -> Property -> a -> ValSet a
setVal (ValSet p v av sl) pos property newVal = 
    let subDict = fromJust $ Map.lookup pos v  
    in ValSet p (Map.insert pos (Map.insert property newVal subDict) v) av sl
    
cartProd:: [[a]] -> [[a]] -> [[a]]
cartProd xs ys = [ x ++ y | x <- xs, y <- ys]

inflowPositions::[Position]
inflowPositions =  makePositions $!  0 :  (map maxPos  $! enumFrom Y) 

makeAllPositions::[Position]
makeAllPositions = makePositions $! map maxPos $! enumFrom X 
         
makePositions :: [Int] -> [Position]
makePositions maxes =
    let ranges = map (\x -> map (:[]) [0..x] ) maxes 
        posCoords = foldl' cartProd [[]] ranges
    in map (\coords -> Position (head coords) (coords!!1) (coords!!2) 0) posCoords
              
--mergeValSets :: (Num a, Fractional a) => ValSet a -> ValSet a-> ValSet a
mergeValSets modifying base = foldl'
    (\prev next-> 
        foldl'
            (\p1 n1-> setVal p1 n1 next $! prop next n1 Center modifying) 
            prev
            $! calculatedPositions modifying 
    )
    base
    $! enumFrom U

modifyPositionComponent::Position -> Direction -> Int -> Position
modifyPositionComponent (Position x y z t) direction amt= case direction of 
    X -> Position amt y z t
    Y -> Position x amt z t 
    Z -> Position x y amt t
    Time -> Position x y z amt

isUpperSide:: Side -> Bool
isUpperSide side = case side of
    East -> True
    North -> True
    Top -> True
    Now -> True 
    _->False

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
    Center -> undefined 
    
offsetPosition:: Position->Side ->Position
offsetPosition (Position x y z t) side = case side of
    Center -> Position x y z t 
    Now -> Position x y z t 
    Prev -> Position x y z $! mod (t - 1) storedSteps  
    _ -> 
        let position = Position x y z t
            maxOrMin = if isUpperSide side then min else max
            offsetAmount = if isUpperSide side then 1 else (-1)
            direction = directionFromCenter side
            boundary = if isUpperSide side 
                then maxPos direction
                else 0
        in modifyPositionComponent position direction 
            $ maxOrMin boundary $ getPositionComponent position direction + offsetAmount   

getPositionComponent:: Position -> Direction -> Int
getPositionComponent (Position x y z t) d = case d of 
    X -> x
    Y -> y
    Z -> z
    Time -> t

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = fromIntegral $! length terms 
        f p n= p + n / len
    in foldl' f 0 terms

positionIfWall (Position x y z t) = if Set.member (Position x y z 0) wallPositionsSet
    then Position x y z 0
    else Position x y z t
    
envIfWall (Position x y z _) env = if Set.member (Position x y z 0) wallPositionsSet
    then id $! initialGrid
    else env       

--prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
prop property position side env = 
    let neighbor = offsetPosition position side
        noValError = error ("no value "
                            ++ show (xPos position)++ " "
                            ++ show (yPos position)++ " "
                            ++ show (zPos position)++ " "
                            ++ show (timePos position)++ " "
                            ++ show property ++ " "
                            ++ show side)
        getVal:: Position -> Map.Map Position (Map.Map Property Double) -> Double
        getVal p set = fromMaybe 
            --(case timePos position of
             --   0 -> noValError
              --  _ -> prop property (offsetPosition p Prev) side  env)
            --noValError
            (case property of
                Density -> noValError
                Temperature -> noValError
                Pressure -> noValError 
                _ -> fromJust $! Map.lookup (modifyPositionComponent p Time 0) (vals initialGrid )>>= Map.lookup property
            )
            --(fromJust $! Map.lookup (offsetPosition p Prev) set >>= Map.lookup property)   
            (Map.lookup p set >>= Map.lookup property)
        --res p = getVal (positionIfWall p) (vals $! envIfWall p env )
        res p = getVal p (vals env )
    in average [ res position, res neighbor]

