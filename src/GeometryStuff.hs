module GeometryStuff where

import qualified Data.Set as Set
import Data.List
import SolutionDomain 
import CalculusUtils

type FlowOb = ([Position],Position)

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
        
connectTwo p1 p2 allPos = 
    let isBetween testPt = 
            let aSq = distanceSquared p1 testPt
                bSq = distanceSquared p2 testPt
                cSq = distanceSquared p1 p2
                x = sqrt $ (aSq + bSq - cSq ) / 2
                cornerDist = sqrt $ 
                    foldl' (\prev next ->
                            let sl = sideLength next testPt 
                            in prev + (sl*sl) 
                        ) 1 (enumFrom X)
            in cornerDist >= x
    in filter isBetween allPos

connectBounds :: [Position] -> [Position] -> [Position]
connectBounds allPos points = fst $ foldl' 
    (\(prevPts, lastPt) next -> ( prevPts ++ connectTwo lastPt next allPos, next) ) 
    ([],last points) points

-- for 3D shapes (and maybe more) you would give the corner points then they would be connected and then all the points on those new edges are connected to each other etc... 