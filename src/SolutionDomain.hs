{-# LANGUAGE ScopedTypeVariables, RankNTypes #-}
module SolutionDomain where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Control.Monad
import Data.Maybe
import Data.List
import GeometryStuff
import System.Random.Shuffle
import System.Random

data Sign = Positive | Zero | Negative deriving (Enum,Ord,Eq,Show)

data Side =  Now | Prev |East | West | North | South | Top | Bottom | Center deriving (Show,Eq,Ord, Enum)
  
data SchemeType = Directional | Nondirectional       
  
data Direction = Time | X | Y | Z deriving (Enum,Ord,Eq,Show)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    lhs:: ![a]
    ,rhs:: ![a]
    ,unknownProperty:: !Property}    
data IntegralType = Body | Surface deriving (Eq)
data Property = Speed | Vorticity | U | V | W | Density | Temperature | Mew | Pressure deriving (Ord,Eq,Enum,Show)
data Position = Position {
    spatialPos:: ![Int]
    ,spatialDimens :: !Int 
    ,timePos:: !Int } deriving (Eq, Ord, Show)
data ValSet a = ValSet{
    calculatedPositions:: ![Position]
    ,vals:: !(Map.Map Position (Map.Map Property a))
    ,areaVal:: !(Map.Map Position (Map.Map Side a))
    ,sideLen:: !(Map.Map Position (Map.Map Direction a))
    ,timeLevel :: Int }
data Expression a = Expression{getTerms:: ![Term a]} 
data Term a = 
    Constant {val:: !a}
    | Unknown { coeff:: !a }
    | SubExpression {expression:: !(Expression a) }  
    | Derivative { denom:: !Direction ,function:: !(SchemeType -> Position->Side->Term a), centered:: !Side, 
        multiplier:: !(SchemeType->Position->Side-> a) } 
        
instance Functor Term  where
    fmap f x = case x of
         Constant c -> Constant $ f c
         Unknown u -> Unknown $ f u
         _ -> error "can't fmap other than unknown and constant"

getSign d = case d of   
    0 -> Zero
    _ -> if d > 0 then Positive else Negative

isConvected :: Property -> Bool
isConvected p = case p of
    Density -> False
    Mew -> False
    Pressure -> False
    --Temperature -> False
    _ -> True

isMomentum :: Property -> Bool
isMomentum p = elem p $ enumFrom U \\ enumFrom Density

convectionFromDirection :: Direction -> Property 
convectionFromDirection d = case d of
    X -> U
    Y -> V
    Z -> W
    _ -> error "no convection for that direction"    
 
directionFromConvection d = case d of
    U -> X
    V -> Y
    W -> Z
    _ -> error "no convection for that direction"    

specificHeatCv :: Double
specificHeatCv = 716

storedSteps:: Int
storedSteps = 3

maxPos :: Direction -> Int
maxPos  d = case d of 
    X -> 180
    Y -> 90
    Z -> 0
    Time -> error "no max time position"
    
gridSize :: Direction -> Double
gridSize d = case d of 
    X -> 1.6
    Y -> 0.8
    Z -> 1
    Time -> error "gridsize for time is the timestep"

cellLength d = gridSize d / (fromIntegral $ maxPos d + 1)

getDirectionComponentIndex direction = case direction of
    X -> 0
    Y -> 1
    Z -> 2
    _ -> error "not defined or not added dimension"

coordToPos coords time = 
    let dirs = enumFrom X
    in Position 
        (map (\x-> round $ coords!!x / (cellLength $ dirs!!x)  ) $ [0.. length coords -1] )
        (length dirs)
        time

posToCoord (Position s _ _ ) direction = 
    cellLength direction * (fromIntegral $ s!!(getDirectionComponentIndex direction)) 

boundaryPair:: Direction -> (Side,Side)
boundaryPair d = case d of 
     X -> (East,West)
     Y -> (North,South)
     Z -> (Top,Bottom)
     Time -> (Now,Now) -- we don't use the previous time step, actually it's more like (future ,now)  but we don't have a future

faceToDirections :: Side -> [Direction] 
faceToDirections s = case s of
    East -> [Y,Z]
    West -> [Y,Z]
    North -> [X,Z]
    South -> [X,Z]
    Top -> [X,Y]
    Bottom -> [X,Y]
    Now -> error "no directions outline the now face"
    Prev -> error "no directions outline the prev face"
    Center -> enumFrom X

removeItems  :: (Ord a, Eq a)=> [a] -> [a]-> [a]
removeItems orig remove= 
    let removeSet = Set.fromList remove
    in filter ( `Set.notMember` removeSet) orig      

boundaryPositions :: [Position]
boundaryPositions = obstacles ++ inflowPositions

boundaryPositionsSet :: Set.Set Position
boundaryPositionsSet = Set.fromList $! boundaryPositions 

obstaclesSet :: Set.Set Position
obstaclesSet = Set.fromList $! obstacles  

obstacle :: Position
obstacle = coordToPos [gridSize X / 4 , gridSize Y / 2 , 0 ] 0

obstacle2 = Position [(quot (maxPos X)  4 ),  (quot (maxPos Y) 2 )+ 5, 0] 3 0

obstacle3 = Position [(quot (maxPos X)  4 ),  (quot (maxPos Y) 2 )+ 2, 0] 3 0

squareBoundsPts :: [Position]
squareBoundsPts = [
   -- obstacle,
   -- offsetPosition (coordToPos [gridSize X / 4 , gridSize Y / 2 , 0 ] 0) West,
    --coordToPos [gridSize X / 4 - 0.12, gridSize Y / 2  , 0 ] 0
     coordToPos [gridSize X / 4 , gridSize Y / 2 + 0.05 , 0 ] 0
    , coordToPos [gridSize X / 4 + 0.1, gridSize Y / 2 + 0.05 , 0 ] 0
    --, coordToPos [gridSize X / 4 + 0.12, gridSize Y / 2  , 0 ] 0
    , coordToPos [gridSize X / 4 + 0.1, gridSize Y / 2 - 0.05 , 0 ] 0
    , coordToPos [gridSize X / 4 , gridSize Y / 2 - 0.05 , 0 ] 0
    ] 

squareBounds :: [Position] 
squareBounds = connectBounds makeAllPositions squareBoundsPts

obstacles :: [Position]
obstacles = 
    let filled =  fillObInterior (Set.fromList squareBounds) $! [offsetPosition obstacle East]
        filledGaps = fillEnclosedGaps makeAllPositions $ Set.fromList filled
    in  -- squareBoundsPts
        --squareBounds
        filled ++ filledGaps

timeStep :: Double            
timeStep = 0.00001 --0.0000001

initialMew =  0.1-- 0.000018

initialTemperature = 223

sutherlandConstant = 120

sutherlandLambda = 
    initialMew *( initialTemperature + sutherlandConstant) / (initialTemperature**1.5)

gasConstantR :: Double
gasConstantR = specificHeatCp - specificHeatCv

specificHeatCp :: Double
specificHeatCp = 1005

initialGridPre:: ValSet Double
initialGridPre= 
    let vMap = foldl' (\prev next -> Map.insert next 
            (case next of 
                U-> 100
                V-> 0
                W-> 0
                Density -> 0.3 -- 1.2
                Pressure -> 10000 -- 101325   doesn't matter, gets solved based on density and temperature. 
                Mew -> initialMew
                Temperature -> initialTemperature
            ) 
            prev) Map.empty (enumFrom U)
        avMap = foldl' (\prev next ->
                    Map.insert next (
                        foldl' (\prev2 next2 -> prev2 * (fromJust $ Map.lookup next2 slMap) ) 1 $ faceToDirections next 
                    ) $! prev
                ) Map.empty $!  enumFrom East
        slMap = foldl' (\prev next -> 
                    Map.insert next (cellLength next ) $! prev
                ) Map.empty $! enumFrom X
        v = foldl' (\prev next -> Map.insert next vMap $! prev) Map.empty  $!  makeAllPositions
        av = foldl' (\prev next -> Map.insert next avMap $! prev) Map.empty $! makeAllPositions
        sl = foldl' (\prev next -> Map.insert next slMap $! prev) Map.empty $! makeAllPositions
    in ValSet makeAllPositions v av sl 0 
    
initialGrid = 
    let calcPos = removeItems (calculatedPositions initialGridPre) $! boundaryPositions
        v= vals initialGridPre
        av = areaVal initialGridPre
        sl = sideLen initialGridPre
    in foldl'
        (\prev nextPos ->
            foldl'
                (\prev2 nextProp -> setVal prev2 nextPos nextProp 0.0) 
                prev
                $ enumFrom U \\ enumFrom Density     
        )
        (ValSet calcPos v av sl 0)
         $! obstacles
    
setVal:: ValSet Double -> Position -> Property -> Double -> ValSet Double
setVal (ValSet p v av sl tl) pos property newVal = 
    let subDict = case Map.lookup pos v  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p (Map.insert pos (Map.insert property newVal subDict) v) av sl tl

setFaceArea :: ValSet Double -> Position -> Side -> Double -> ValSet Double
setFaceArea (ValSet p v av sl tl) pos side newVal = 
    let subDict = case Map.lookup pos av  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v (Map.insert pos (Map.insert side newVal subDict) av) sl tl             
    
setCellLength :: ValSet Double -> Position -> Direction -> Double -> ValSet Double
setCellLength (ValSet p v av sl tl) pos dir newVal = 
    let subDict = case Map.lookup pos sl  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v av (Map.insert pos (Map.insert dir newVal subDict) sl) tl
        
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

setElem newElem index list = 
    map (\x -> if x == index then newElem else list!!x) [0..length list -1]

modifyPositionComponent::Position -> Direction -> Int -> Position
modifyPositionComponent (Position p d t) direction amt= case direction of 
    Time -> Position p d amt
    _ -> Position (setElem amt (getDirectionComponentIndex direction) p) d t

offsetPositionComponent pos dir amt =
    let curVal = getPositionComponent pos dir
    in modifyPositionComponent pos dir (curVal+amt)

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
    Center -> error "there is no direction from center" 
    
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

getPositionComponent:: Position -> Direction -> Int
getPositionComponent (Position p d t) direction = case direction of 
    Time -> t
    _ -> p!!(getDirectionComponentIndex direction)

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = fromIntegral $! length terms 
        f p n= p + n / len
    in foldl' f 0 terms

isBoundaryPosition (Position p d t) = Set.member (Position p d 0) boundaryPositionsSet

isObstaclePosition (Position p d t) = Set.member (Position p d 0) obstaclesSet 

positionIfWall (Position p d t) = if isBoundaryPosition (Position p d t) 
    then Position p d 0
    else Position p d t
    
envIfWall (Position p d _) env = if isBoundaryPosition (Position p d 0) 
    then id $! initialGrid
    else env       

pecletNumber position side env = 
    let direc = directionFromCenter side
        momentum = propCentralDiff (convectionFromDirection direc) position side env  
        density = propCentralDiff Density position side env
        viscosity = propCentralDiff Mew position side env
        l = sideLength direc position env 
    in (density * momentum * l) / viscosity 

--approximateDerivative::(Num a, Fractional a)=>  Term a -> Position-> ValSet a -> Term a
approximateDerivative deriv position vs= case deriv of 
    (Derivative direction func side m) ->
        if directionFromCenter side == direction 
        then
            let neighbor = offsetPosition position side
                sl =  sideLength direction position vs
                sln = sideLength direction neighbor vs
                interval = average [sl, sln ]
                thisVal = func Nondirectional position Center 
                neighborVal = func Nondirectional neighbor Center
                neighborIsUpper = isUpperSide side  
                f first = case (first, neighborIsUpper) of
                    (True, True) -> neighborVal
                    (False, False) -> neighborVal
                    _-> thisVal
                mult = m Nondirectional position Center
            in case (f True, f False) of
                (Constant c1 , Constant c2) -> Constant $ (c1-c2)*mult/interval 
                _ -> error "can't approx >1 order derivs. deriv must produce constants" 
        else
            let (s1, s2) = boundaryPair direction
                n1 = (offsetPosition position s1)
                n2 = (offsetPosition position s2)
                sl =  sideLength direction position vs
                sln1 = sideLength direction n1 vs
                sln2 = sideLength direction n2 vs
                interval = 2 * average [sl,sln1,sln2]
                n1Val = func Nondirectional n1 side
                n2Val = func Nondirectional n2 side
                mult = m Nondirectional position side
            in case (n1Val, n2Val) of
                (Constant c1 , Constant c2) -> 
                    Constant $ (c1-c2)*mult/interval  
                _ -> error "can't approx >1 order derivs. deriv must produce constants"
    _ -> error "can't approx something that isn't a deriv"

magnitude v = sqrt $ foldl' (\prev nVal ->prev + (nVal*nVal)) 0.0 v

prop schemeType =
    let f scheme = (\property pos side env ->
            if side == Center || side == Now || side == Prev
            then propCentralDiff property pos side env 
            else scheme property pos side env) 
    in 
        let scheme = case schemeType of
                Directional -> f propLimitedSlope
                Nondirectional -> f propLimitedSlope
            momentums = enumFrom U \\ enumFrom Density
        in \property pos side env-> case property of
                Speed-> magnitude $ map (\next -> scheme next pos side env) momentums   
                Vorticity -> 
                    let pairs = [(x,y) | x <- momentums , y <- tail $ dropWhile (/= x) momentums ]
                        vortComponents = map (\(a,b) ->
                                let makeDeriv m n 
                                        = Derivative (directionFromConvection m) 
                                            (\_ p s -> Constant $ scheme n p s env) East -- east is arbitrary here  
                                            (\_ _ _ -> 1) 
                                    deriv1 = makeDeriv a b
                                    deriv2 = makeDeriv b a
                                in (val $ approximateDerivative deriv1 pos env)
                                    - (val $ approximateDerivative deriv2 pos env)    
                            ) pairs
                    in magnitude vortComponents 
                Mew -> 
                    let t = scheme Temperature pos side env 
                    in sutherlandLambda * (t**1.5)/ (t + sutherlandConstant)
                Pressure ->
                    gasConstantR * (scheme Temperature pos side env)
                    * (scheme Density pos side env)
                _-> scheme property pos side env  

propDirectional property position side env =
    let neighbor = offsetPosition position side
        decide =
            let peclet = pecletNumber position side env  
            in if peclet == 0
                then propCentralDiff
                else 
                    if abs peclet < 0.99
                    then propCentralDiff
                    else if abs peclet > 2.6
                        then propUpwindDiff peclet  
                        else propQUICK peclet 
    in case (isObstaclePosition neighbor || isObstaclePosition position
                , isMomentum property 
                ,elem side ( enumFrom East \\ enumFrom Center )
                    && isConvected property ) of
        --(True,True,_)-> 0.0
        (_,_,True) -> decide property position side env 
        _ -> propCentralDiff property position side env 
        
oppositeSide :: Side -> Side 
oppositeSide s = 
    let d = directionFromCenter s
        (s1,s2) = boundaryPair d
    in if isUpperSide s 
        then s2 else s1 

propQUICK :: Double -> Property -> Position -> Side ->ValSet Double -> Double
propQUICK pec property position side env = 
    let upper = if isUpperSide side then side else Center 
        lower = if isUpperSide side then Center else side
        doubleOffset pos = offsetPosition pos >>= offsetPosition
        farPoint upBracket downBracket = if upBracket == Center
            then offsetPosition position $ oppositeSide downBracket
            else doubleOffset position upBracket 
    in if pec >= 0
        then ((6/8)* propCentralDiff property (offsetPosition position lower) Center env ) 
            + ((3/8)* propCentralDiff property (offsetPosition position upper) Center env ) 
            - ((1/8)* propCentralDiff property (farPoint lower upper) Center env ) 
        else ((6/8)* propCentralDiff property (offsetPosition position upper) Center env ) 
            + ((3/8)* propCentralDiff property (offsetPosition position lower) Center env ) 
            - ((1/8)* propCentralDiff property (farPoint upper lower) Center env )

propUpwindDiff :: Double -> Property -> Position -> Side -> ValSet Double -> Double              
propUpwindDiff pec property position side env = 
    let upper = if isUpperSide side then side else Center 
        lower = if isUpperSide side then Center else side
    in if pec >= 0
        then propCentralDiff property (offsetPosition position lower) Center env
        else propCentralDiff property (offsetPosition position upper) Center env

-- http://www.ammar-hakim.org/_static/files/1010-muscl-hancock.pdf
-- http://bulletin.pan.pl/(60-1)45.pdf
propLimitedSlope property position side env = 
    let valCenter = propCentralDiff property position Center env
        d = directionFromCenter side
        interval = sideLength d position env
        (upper, lower) = boundaryPair d
        upperNVal:(lowerNVal:_)  
            = map (\s -> propCentralDiff property (offsetPosition position s) Center env) 
                [upper,lower] 
        --ave = superbee (upperNVal - valCenter) (valCenter - lowerNVal)
        --ave = epsilon (upperNVal - valCenter) (valCenter - lowerNVal) (interval *interval *interval)
        --ave = minmodLimit (upperNVal - valCenter) (valCenter - lowerNVal)
        ave = ospre (upperNVal - valCenter) (valCenter - lowerNVal)
        -- ave = vanLeer (upperNVal - valCenter) (valCenter - lowerNVal)
        --ave = ((upperNVal - valCenter)+(valCenter - lowerNVal))/2
    in (if isUpperSide side then (+) else (-)) valCenter  
            $ 0.5 * ave   

--prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
propCentralDiff property position side env = 
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
            --noValError
            (case property of
                Density -> noValError
                Temperature -> noValError
                Pressure -> noValError 
                _ -> case  Map.lookup (modifyPositionComponent p Time 0) (vals initialGrid )>>= Map.lookup property of
                        Nothing -> noValError
                        Just r -> r
            )   
            (Map.lookup p set >>= Map.lookup property)
        res p = getVal p (vals env )
    in case (isObstaclePosition neighbor || isObstaclePosition position
                , isMomentum property, position == neighbor) of
        --(True,True,_) -> 0.0
        (_,_,True) -> res position
        _ -> average [ res position, res neighbor]

-- GEOMETRY STUFF 
 
fillObInterior :: Set.Set Position -> [Position] -> [Position]
fillObInterior obBounds points = case points of
    [] -> Set.toList obBounds
    x:xs -> 
        let unvN = filter (`Set.notMember` obBounds) $ getNeighbors x
        in fillObInterior (Set.union obBounds $ Set.fromList unvN) $ xs ++ unvN
        
distanceSquared :: Position -> Position -> Double        
distanceSquared (Position sp1 _ _) (Position sp2 _ _) = 
    let components = zip sp1 sp2         
    in fromIntegral $ foldl' (\prev (a,b) -> prev + (a-b)*(a-b) ) 0 components
    
distance:: Position -> Position -> Double    
distance p1 p2 = sqrt $ distanceSquared p1 p2 
        
isBetween p1 p2 testPt vs= 
    let aSq = distanceSquared p1 testPt
        bSq = distanceSquared p2 testPt
        cSq = distanceSquared p1 p2
        x = sqrt $ (aSq + bSq - cSq ) / 2
        cornerDist = sqrt $ 
            foldl' (\prev next ->
                    let sl = sideLength next testPt vs
                    in prev + (sl*sl) 
                ) 0 (enumFrom X)
    in cornerDist > x
        && testPt /= p1 && testPt /= p2        
        
connectTwo p1 p2 allPos vs = filter (isBetween p1 p2 vs) allPos

getNeighborsCorners position =
    let cross = filter (\(a,b) -> a /= b ) [(x,y) | x <- enumFrom East, y<-enumFrom East] 
        n = nub $ filter (\x -> x /= position)
            $ map (\(a,b) -> offsetPosition (offsetPosition position a) b ) cross
    in shuffle' n (length n) $ mkStdGen $ sum $ spatialPos position

getNeighbors position =
    nub $ filter (\x -> x /= position)
        $ map (\x -> offsetPosition position x ) $ enumFrom East

tracePath g end path = 
    let prev = Map.lookup end g
    in case prev of
        Nothing -> path
        Just (_, p) -> if p == end 
            then path 
            else tracePath g p $ p:path  

hugeNumber = 99999999999.9
rounder f n = (fromInteger $ round $ f * (10^n)) / (10.0^^n)

shortestPath unvisited end visited = if unvisited==Map.empty 
    then error "problem in shortest path.. searched everything and wasn't found" 
    else 
        let comparer a= if mod (Map.size unvisited * Map.size visited + a ) 2 < 1 
                then (<) else (<=)
            (pos, (dist,predecessor )) = 
                foldl' (\(p1,(p2,p3) ) (n1,(n2,n3)) -> if comparer (fromIntegral $ round $ n2*p2*57) n2  p2 
                            then (n1,(n2,n3)) else (p1,(p2,p3))  ) 
                        (end,(hugeNumber,end)) 
                        $ shuffle' (Map.toList unvisited) (Map.size unvisited) $ mkStdGen 137
            neighbors = filter (\y -> Map.notMember y visited) $ getNeighborsCorners pos 
            updatedUnvisited = foldl' (\prev next -> 
                    let newd = (distance next pos) + dist
                        oldd = fst $ fromMaybe (hugeNumber,end) (Map.lookup next prev) 
                    in if comparer (fromIntegral $ round $ newd*oldd*57) newd oldd 
                        then Map.insert next (newd,pos) prev
                        else prev
                ) (Map.delete pos unvisited) neighbors
            updatedVisited = Map.insert pos (dist,predecessor ) visited
        in case pos == end of
            True -> tracePath updatedVisited end []
            False -> shortestPath updatedUnvisited end updatedVisited
           
connectBounds allPos points = fst $ foldl' 
    (\(prevPts, lastPt) next -> 
        ( prevPts 
            ++ shortestPath 
                (Map.insert lastPt (0.0,lastPt) Map.empty) 
                next 
                Map.empty 
        , next) ) 
    ([],last points) points

fillEnclosedGaps :: [Position] -> Set.Set Position -> [Position]
fillEnclosedGaps allPos wallPos = 
    let added = filter
            (\x ->
                let opposingWalls = foldl'
                        (\prev next ->
                            let (s1 ,s2) = boundaryPair next
                            in prev ||
                                ( Set.member (offsetPosition x s1) wallPos 
                                    && Set.member (offsetPosition x s2) wallPos )
                        )
                        False
                        $ enumFrom X
                in (Set.notMember x wallPos)
                    && opposingWalls
            )
            allPos
    in case added of 
        [] -> []
        _ -> added ++ (fillEnclosedGaps allPos $ Set.union wallPos $ Set.fromList added) 

-- sideArea:: (Num a, Fractional a)=>Side -> Position -> a
sideArea s (Position p d _) vs = case s of 
    Now -> 1
    Prev -> 1
    _ -> case  Map.lookup (Position p d 0) (areaVal $! vs)  >>= Map.lookup s of
            Nothing -> error $ "error getting side area for " ++ show p ++ " " ++ show s
            Just sa -> sa

-- sideLength:: (Num a, Fractional a) => Direction -> Position ->  a
sideLength direction (Position p d _) vs = case direction of 
    Time -> timeStep
    _ -> case   Map.lookup (Position p d 0) (sideLen $! vs) >>= Map.lookup direction of
            Nothing -> error $ "error getting side length for "++ show p++ " " ++ show direction
            Just len -> len    

chooseSlopeHelper f1 f2 x y =
    let sign = getSign x 
    in if (getSign x) == (getSign y)
        then case sign of
                Positive -> Just $ f1 x y
                Negative -> Just $ f2 x y
                Zero -> Nothing 
        else Nothing 

chooseSlope f1 f2 n = 
    let res = foldM
            (\prev next -> chooseSlopeHelper f1 f2 prev next) 
            (n!!0) n
    in case res of
        Just x -> x
        Nothing -> 0.0 

minmod = chooseSlope min max
maxmod = chooseSlope max min

superbee a b = minmod [maxmod [a,b], minmod [2*a,2*b] ]

minmodLimit a b = minmod [ (a + b) /2 , 2*a, 2*b ]

epsilon a b eSq = ( ((b*b + eSq )*a)  + ((a*a+eSq)*b) ) / ((a*a) + (b*b) + (2*eSq)) 

vanLeer = altFormLimiter (\r -> let abs_r = abs r in (r+abs_r)/(1+abs_r) )  

ospre = altFormLimiter (\r -> 1.5 * (r*r + r) / (r*r + r + 1) )

altFormLimiter rFunc a b = 
    let expression steep shallow = 
            let r = steep/shallow
            in shallow * rFunc r
    in case (getSign a, getSign b) of
        (Positive,Positive) -> expression (max a b) (min a b)
        (Negative,Negative) -> expression  (min a b) (max a b)
        _-> 0.0 


       