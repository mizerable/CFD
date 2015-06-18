
module SolutionDomain where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Data.Maybe
import Control.Monad.Reader as Reader
import Data.List
import GeometryStuff

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
    spatialPos:: ![Int]
    ,spatialDimens :: !Int 
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

timeStep :: Double            
timeStep = 0.005

specificHeatCv :: Double
specificHeatCv = 15

storedSteps:: Int
storedSteps = 4

maxPos:: Direction -> Int
maxPos d = case d of 
    X -> 450
    Y -> 150
    Z -> 0
    Time -> undefined

removeItems  :: (Ord a, Eq a)=> [a] -> [a]-> [a]
removeItems orig remove= 
    let removeSet = Set.fromList remove
    in filter ( `Set.notMember` removeSet) orig      

wallPositions :: [Position]
wallPositions = (obstacles ++ inflowPositions)

wallPositionsSet :: Set.Set Position
wallPositionsSet = Set.fromList $! wallPositions 

obstacle :: Position
obstacle = Position [quot (maxPos X)  3,  quot (maxPos Y) 2, 0] 3 0

obstacle2 = Position [(quot (maxPos X)  3 ),  (quot (maxPos Y) 2 )+ 5, 0] 3 0

obstacle3 = Position [(quot (maxPos X)  3 ),  (quot (maxPos Y) 2 )+ 2, 0] 3 0

squareBoundsPts :: [Position]
squareBoundsPts = [
    obstacle
    , Position [(quot (maxPos X)  3 ),  (quot (maxPos Y) 2 )+ 3, 0] 3 0
    , Position [(quot (maxPos X)  3 )+ 3,  (quot (maxPos Y) 2 ) + 3, 0] 3 0
    , Position [(quot (maxPos X)  3 ) + 6 ,  (quot (maxPos Y) 2 )+ 3, 0] 3 0
    , Position [(quot (maxPos X)  3) + 6,  (quot (maxPos Y) 2), 0] 3 0
    , Position [(quot (maxPos X)  3) + 6,  (quot (maxPos Y) 2 )- 3, 0] 3 0
    , Position [(quot (maxPos X)  3) + 3,  (quot (maxPos Y) 2) - 3, 0] 3 0
    , Position [(quot (maxPos X)  3),  (quot (maxPos Y) 2 - 3), 0] 3 0
    ] 

squareBounds :: [Position] 
squareBounds = connectBounds makeAllPositions squareBoundsPts

obstacles :: [Position]
obstacles = squareBounds ++ (fillObInterior (Set.fromList squareBounds) $! offsetPosition obstacle East)

initialGridPre:: ValSet Double
initialGridPre= 
    let vMap = foldl' (\prev next -> Map.insert next 
            (case next of 
                U-> 2.25
                V-> 0
                W-> 0
                Density -> 1
                _-> 1
                ) 
            prev) Map.empty (enumFrom U)
        avMap = foldl' (\prev next -> Map.insert next 1 $! prev) Map.empty $!  enumFrom East
        slMap = foldl' (\prev next -> Map.insert next 1 $! prev) Map.empty $! enumFrom X
        v = foldl' (\prev next -> Map.insert next vMap $! prev) Map.empty  $!  makeAllPositions
        av = foldl' (\prev next -> Map.insert next avMap $! prev) Map.empty $! makeAllPositions
        sl = foldl' (\prev next -> Map.insert next slMap $! prev) Map.empty $! makeAllPositions
    in ValSet makeAllPositions v av sl 
    
initialGrid = 
    let calcPos = removeItems (calculatedPositions initialGridPre) $! wallPositions
        v= vals initialGridPre
        av = areaVal initialGridPre
        sl = sideLen initialGridPre
    in foldl'
        (\prev next -> setVal prev next U 0.0)
        (ValSet calcPos v av sl)
         $! obstacles
        
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
    in map (\coords -> Position coords (length coords) 0) posCoords
              
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

setElem newElem index list = 
    map (\x -> if x == index then newElem else list!!x) [0..length list -1]

modifyPositionComponent::Position -> Direction -> Int -> Position
modifyPositionComponent (Position p d t) direction amt= case direction of 
    Time -> Position p d amt
    _ -> Position (setElem amt (getDirectionComponentIndex direction) p) d t

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
offsetPosition (Position p d t) side = case side of
    Center -> Position p d t 
    Now -> Position p d t 
    Prev -> Position p d $! mod (t - 1) storedSteps  
    _ -> 
        let position = Position p d t
            maxOrMin = if isUpperSide side then min else max
            offsetAmount = if isUpperSide side then 1 else (-1)
            direction = directionFromCenter side
            boundary = if isUpperSide side 
                then maxPos direction
                else 0
        in modifyPositionComponent position direction 
            $ maxOrMin boundary $ getPositionComponent position direction + offsetAmount   

getDirectionComponentIndex direction = case direction of
    X -> 0
    Y -> 1
    Z -> 2
    _ -> error "not defined or not added dimension"

getPositionComponent:: Position -> Direction -> Int
getPositionComponent (Position p d t) direction = case direction of 
    Time -> t
    _ -> p!!(getDirectionComponentIndex direction)

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = fromIntegral $! length terms 
        f p n= p + n / len
    in foldl' f 0 terms

positionIfWall (Position p d t) = if Set.member (Position p d 0) wallPositionsSet
    then Position p d 0
    else Position p d t
    
envIfWall (Position p d _) env = if Set.member (Position p d 0) wallPositionsSet
    then id $! initialGrid
    else env       

--prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
prop property position side env = 
    let neighbor = offsetPosition position side
        noValError = error ("no value "
                            ++ 
                            (
                                foldl' (\prev next -> prev ++ " " ++  (show $ (spatialPos position)!!next )) "" 
                                    [0..spatialDimens position-1] 
                            )++ " "
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

-- GEOMETRY STUFF 

--the bounds are NOT in the result. 
fillObInterior :: Set.Set Position -> Position -> [Position]
fillObInterior obBounds point = 
    let neighbors = map (offsetPosition point) $ enumFrom East
        unvisitedNeighbors = 
            filter (\x-> x /= point && Set.notMember x obBounds) neighbors
    in point --i think the strictness of foldl' is important here
        : foldl' (\prev next -> prev ++ fillObInterior (Set.insert point obBounds) next) 
            [] unvisitedNeighbors
        
distanceSquared :: Position -> Position -> Double        
distanceSquared p1 p2 = 
    fromIntegral $
    foldl' (\prev next -> 
            let diff = (getPositionComponent p1 next) - (getPositionComponent p2 next) 
            in prev + ( diff * diff )      
        ) 0 $ enumFrom X         
    
distance p1 p2 = sqrt $ distanceSquared p1 p2 
        
isBetween p1 p2 testPt = 
    let aSq = distanceSquared p1 testPt
        bSq = distanceSquared p2 testPt
        cSq = distanceSquared p1 p2
        x = sqrt $ (aSq + bSq - cSq ) / 2
        cornerDist = sqrt $ 
            foldl' (\prev next ->
                    let sl = sideLength next testPt 
                    in prev + (sl*sl) 
                ) 0 (enumFrom X)
    in cornerDist >= x        
        
connectTwo p1 p2 allPos = filter (isBetween p1 p2) allPos

connectBounds :: [Position] -> [Position] -> [Position]
connectBounds allPos points = fst $ foldl' 
    (\(prevPts, lastPt) next -> ( prevPts ++ connectTwo lastPt next allPos, next) ) 
    ([],last points) points


-- sideArea:: (Num a, Fractional a)=>Side -> Position -> a
sideArea s (Position p d _) = case s of 
    Now -> 1
    Prev -> 1
    _ -> fromJust $! Map.lookup (Position p d 0) (areaVal $! initialGridPre)  >>= Map.lookup s    

-- sideLength:: (Num a, Fractional a) => Direction -> Position ->  a
sideLength direction (Position p d _) = case direction of 
    Time -> timeStep
    _ -> fromJust $! Map.lookup (Position p d 0) (sideLen $! initialGridPre) >>= Map.lookup direction   

