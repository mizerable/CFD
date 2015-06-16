module GeometryStuff where

type TwoDShape = [(Int,Int)]
type TwoDPoint = (Int,Int)

data Quadrant = NE | SE | NW | SW 

isInterior :: TwoDShape -> TwoDPoint -> Bool
isInterior shapeList (x,y) = undefined


getSector :: TwoDPoint -> TwoDPoint -> TwoDPoint -> Int 
getSector (x1,y1) (x2, y2) (xp, yp) = 
    let westPoint = if x1 < x2 then (x1,y1) else (x2,y2)
        eastPoint = if fst westPoint == x1 then (x2,y2) else (x1,y1) 
    in undefined
    
combineQuadrants :: Quadrant  -> Quadrant -> Int
combineQuadrants eastQ westQ = undefined