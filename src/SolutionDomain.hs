{-# LANGUAGE ScopedTypeVariables, RankNTypes #-}
module SolutionDomain where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Data.Maybe
import Data.List

data Side =  Now | Prev |East | West | North | South | Top | Bottom | Center deriving (Show,Eq,Ord, Enum)
  
data SchemeType = Directional | Nondirectional       
  
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
    | Derivative { denom:: !Direction ,function:: !(SchemeType -> Position->Side->Term a), centered:: !Side, 
        multiplier:: !(SchemeType->Position->Side-> a) } 
        
instance Functor Term  where
    fmap f x = case x of
         Constant c -> Constant $ f c
         Unknown u -> Unknown $ f u
         _ -> error "can't fmap other than unknown and constant"

isConvected :: Property -> Bool
isConvected p = case p of
    --Density -> False
    Mew -> False
    --Pressure -> False
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

timeStep :: Double            
timeStep = 0.01

specificHeatCv :: Double
specificHeatCv = 15

storedSteps:: Int
storedSteps = 4

maxPos :: Direction -> Int
maxPos  d = case d of 
    X -> 300
    Y -> 100
    Z -> 0
    Time -> error "no max time position"
    
gridSize :: Direction -> Double
gridSize d = case d of 
    X -> 300
    Y -> 100
    Z -> 1
    Time -> error "gridsize for time is the timestep"

cellLength d = gridSize d / (fromIntegral $ maxPos d + 1)

coordToPos coords time = 
    let dirs = enumFrom X
    in Position 
        (map (\x-> round $ coords!!x / (cellLength $ dirs!!x)  ) $ [0.. length coords -1] )
        (length dirs)
        time

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
    coordToPos [gridSize X / 4 , gridSize Y / 2 + 5 , 0 ] 0
    , coordToPos [gridSize X / 4 + 10, gridSize Y / 2 + 5 , 0 ] 0
    , coordToPos [gridSize X / 4 + 10, gridSize Y / 2 - 5 , 0 ] 0
    , coordToPos [gridSize X / 4 , gridSize Y / 2 - 5 , 0 ] 0
    ] 

squareBounds :: [Position] 
squareBounds = connectBounds makeAllPositions squareBoundsPts

obstacles :: [Position]
obstacles = 
    let filled =  fillObInterior (Set.fromList squareBounds) $! [offsetPosition obstacle East]
        filledGaps = fillEnclosedGaps makeAllPositions $ Set.fromList filled
    in  --squareBoundsPts
        --squareBounds
        filled ++ filledGaps

initialGridPre:: ValSet Double
initialGridPre= 
    let vMap = foldl' (\prev next -> Map.insert next 
            (case next of 
                U-> 2
                V-> 0
                W-> 0
                Density -> 1
                _-> 1
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
    in ValSet makeAllPositions v av sl 
    
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
        (ValSet calcPos v av sl)
         $! obstacles
         
upScaleGrid :: ValSet Double -> Int ->Direction -> ValSet Double
upScaleGrid smallGrid scale direction = 
    foldl'
        (\vs nextPos ->
            let newp =  upScalePosition nextPos scale direction
                (ValSet p newv newav newsl ) = 
                    (\env -> upScaleVals nextPos smallGrid env newp)
                    .(\env -> upScaleAV nextPos smallGrid env newp direction $ fromIntegral scale)
                    .(\env -> upScaleSL nextPos smallGrid env newp direction $ fromIntegral scale) 
                         $ vs
            in ValSet (newp ++ p) newv newav newsl -- tail because the first is duplicated    
        ) emptyValSet
        (calculatedPositions smallGrid)

upScalePosition :: Position -> Int -> Direction -> [Position]
upScalePosition (Position p d t) scale direction = 
    let idx = getDirectionComponentIndex direction
        anchor = p!!idx
    in map (\x -> Position (setElem (anchor*scale + x) idx p) d t) [0.. scale-1]
    
upScaleVals:: Position -> ValSet Double -> ValSet Double -> [Position] -> ValSet Double   
upScaleVals oldPos oldGrid vs newPoses = foldl'
    (\prev1 nextProp ->
        foldl'
            (\prev2 nextPos -> setVal prev2 nextPos nextProp 
                $ fromMaybe (error $ "couldn't get while upscaling " ++ show oldPos ++" "++ show nextProp)
                    $ Map.lookup oldPos (vals oldGrid) >>= Map.lookup nextProp)
            prev1
            newPoses   
    ) vs $ enumFrom U

upScaleAV :: Position -> ValSet Double -> ValSet Double -> [Position] -> Direction -> Double-> ValSet Double  
upScaleAV oldPos oldGrid vs newPoses direction scale= 
    let (s1, s2) = boundaryPair direction
    in foldl'
        (\prev nextSide -> 
            foldl'
                (\prev2 nextPos ->
                    let changeScale = if nextSide == s1 || nextSide == s2
                            then 1 else scale -- for example, if i increase x resolution, then only east/west faces are unchanged
                    in setFaceArea prev2 nextPos nextSide ( sideArea nextSide oldPos oldGrid / changeScale) 
                ) prev newPoses
        ) vs $ enumFrom East
    
upScaleSL :: Position ->  ValSet Double ->ValSet Double -> [Position] -> Direction -> Double ->ValSet Double
upScaleSL oldPos oldGrid vs newPoses direction scale = 
    foldl'
        (\prev0 nextDir ->
            foldl'
                (\prev nextPos -> 
                    let changeScale = if nextDir == direction then scale else 1
                    in setCellLength prev nextPos nextDir $ sideLength nextDir oldPos oldGrid / changeScale
                ) prev0 $ newPoses
        ) vs $ enumFrom X
    
emptyValSet :: ValSet Double    
emptyValSet = ValSet [] Map.empty Map.empty Map.empty
    
setVal:: ValSet Double -> Position -> Property -> Double -> ValSet Double
setVal (ValSet p v av sl) pos property newVal = 
    let subDict = case Map.lookup pos v  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p (Map.insert pos (Map.insert property newVal subDict) v) av sl

setFaceArea :: ValSet Double -> Position -> Side -> Double -> ValSet Double
setFaceArea (ValSet p v av sl) pos side newVal = 
    let subDict = case Map.lookup pos av  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v (Map.insert pos (Map.insert side newVal subDict) av) sl             
    
setCellLength :: ValSet Double -> Position -> Direction -> Double -> ValSet Double
setCellLength (ValSet p v av sl) pos dir newVal = 
    let subDict = case Map.lookup pos sl  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v av (Map.insert pos (Map.insert dir newVal subDict) sl)
        
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

prop schemeType = case schemeType of
    Directional -> propDirectional
    Nondirectional -> propCentralDiff

propDirectional property position side env =
    let neighbor = offsetPosition position side
        decide =
            let peclet = pecletNumber position side env  
            in if abs peclet > 99999999
                then propUpwindDiff peclet  
                else propQUICK peclet 
    in case (isObstaclePosition neighbor || isObstaclePosition position
                , isMomentum property 
                ,elem side ( enumFrom East \\ enumFrom Center )
                    && isConvected property ) of
        (True,True,_)-> 0.0
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
        (True,True,_) -> 0.0
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
distanceSquared p1 p2 = 
    fromIntegral $
    foldl' (\prev next -> 
            let diff = (getPositionComponent p1 next) - (getPositionComponent p2 next) 
            in prev + ( diff * diff )      
        ) 0 $ enumFrom X         
        
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

getNeighbors position =
    map (\x-> offsetPosition position x)  
    $ filter (\x -> position /= (offsetPosition position x) ) $ enumFrom East

tracePath g end path = 
    let prev = Map.lookup end g
    in case prev of
        Nothing -> path
        Just p -> if p == end 
            then path 
            else tracePath g p $ p:path  

shortestPath :: [Position] -> Position -> Map.Map Position Position -> [Position]
shortestPath q end g = case q of
    [] -> error "problem in shortest path.. searched everything and wasn't found" 
    x:xs ->  case x == end of 
        True -> tracePath g end []
        False ->
            let n = filter (\y -> Map.notMember y g) $ getNeighbors x
            in shortestPath (xs ++ n) end 
                $ foldl' (\prev next -> Map.insert next x prev ) g n  
            
connectBounds :: [Position] -> [Position] -> [Position]
connectBounds allPos points = fst $ foldl' 
    (\(prevPts, lastPt) next -> ( prevPts ++ shortestPath [lastPt] next (Map.insert lastPt lastPt Map.empty) , next) ) 
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

