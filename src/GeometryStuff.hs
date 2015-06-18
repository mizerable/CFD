module GeometryStuff where

import qualified Data.Set as Set
import SolutionDomain 

type FlowOb = ([Position],Position)

--the bounds are NOT in the result. 
fillObInterior :: Set.Set Position -> Position -> [Position]
fillObInterior obBounds point = 
    let neighbors = map (offsetPosition point) enumFrom East
        unvisitedNeighbors = 
            filter (x-> (x != point) && Set.notMember x obBounds) neighbors
    in point --i think the strictness of foldl' is important here
        : foldl' (\prev next -> prev ++ fillObInterior (Set.inert point obBounds) next) 
            [] unvisitedNeighbors
        
connectTwo p1 p2 allPos = 
    let isBetween testPt = undefined -- involves solving for intersection of lines of sides inside the interval of the cell.
    in filter isBetween allPos

connectBounds :: [Position] -> [Position] -> [Position]
connectBounds allPos points = fst $ foldl' 
    (\(prevPts, lastPt) next -> ( prevPts ++ connectTwo lastPt next allPos, next) ) 
    ([],last points) points

-- for 3D shapes (and maybe more) you would give the corner points then they would be connected and then all the points on those new edges are connected to each other etc... 