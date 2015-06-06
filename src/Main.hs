module Main where

import qualified Data.Map as Map 
import Data.Maybe
import Control.Monad.State as State
import Control.Monad.Reader as Reader
  
data Side = East | West | North | South | Top | Bottom | Now | Prev | Center deriving (Show,Eq,Ord, Enum)
data Direction = Time | X | Y | Z deriving (Enum,Ord,Eq,Show)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    lhs::[a]
    ,rhs::[a]
    ,unknownProperty::Property}    
data IntegralType = Body | Surface deriving (Eq)
data Property = U | V | W | Density | Temperature | Mew deriving (Ord,Eq,Enum,Show)
data Position = Position {
    xPos::Int
    ,yPos::Int
    ,zPos::Int
    ,timePos::Double } deriving (Eq, Ord, Show)
data ValSet a = ValSet{
    positions:: [Position]
    ,vals:: Map.Map Position (Map.Map Property a)
    ,areaVal::Map.Map Position (Map.Map Side a)
    ,sideLen:: Map.Map Position (Map.Map Direction a) }
data Expression a = Expression{getTerms::[Term a]} 
data Term a = Constant {val::a} | Unknown { coeff::a } | SubExpression {expression::Expression a} 
    | Derivative { denom::Direction ,function:: Position->Side->Term a, centered::Side }

timeStep::Double
timeStep = 0.0001 

maxPos:: Direction -> Int
maxPos d = case d of 
    X -> 12
    Y -> 12
    Z -> 12

getPositionComponent (Position x y z _) d = case d of 
    X -> x
    Y -> y
    Z -> z

modifyPositionComponent (Position x y z t) direction amt= case direction of 
    X -> Position (amt) y z t
    Y -> Position x (amt) z t 
    Z -> Position x y (amt) t
    
offsetPosition:: Position->Side ->Position
offsetPosition (Position x y z t) side = case side of
    Center -> Position x y z t 
    Now -> Position x y z t 
    Prev -> Position x y z (max 0 (t - timeStep))  
    _ -> 
        let position = Position x y z t
            maxOrMin = if isUpperSide side then min else max
            offsetAmount = if isUpperSide side then 1 else (-1)
            direction = directionFromCenter side
            boundary = if isUpperSide side 
                then maxPos direction 
                else 0
        in modifyPositionComponent position direction 
            $ maxOrMin boundary $ (getPositionComponent position direction) + offsetAmount   
        
prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
prop property position side = do
    (ValSet _ v _ _) <- ask 
    return $ 
        let neighbor = offsetPosition position side
            getVal p = fromJust $ Map.lookup p v >>= Map.lookup property  
        in average [getVal position,getVal neighbor]

initialGrid::(Num a) => ValSet a
initialGrid= 
    let p = makePositions
        vMap = foldr (\next -> \prev -> Map.insert next 1 prev) Map.empty (enumFrom U)
        avMap = foldr (\next -> \prev -> Map.insert next 1 prev) Map.empty (enumFrom East)
        slMap = foldr (\next -> \prev -> Map.insert next 1 prev) Map.empty (enumFrom X)
        v = foldr (\next -> \prev -> Map.insert next vMap prev) Map.empty makePositions
        av = foldr (\next -> \prev -> Map.insert next avMap prev) Map.empty makePositions
        sl = foldr (\next -> \prev -> Map.insert next slMap prev) Map.empty makePositions
    in ValSet p v av sl
    
cartProd:: [[a]] -> [[a]] -> [[a]]
cartProd xs ys = [ x ++ y | x <- xs, y <- ys]

makePositions::[Position]
makePositions = 
    let ranges = (enumFrom X) |> map maxPos |> map (\x -> [1..x] |> map (\y -> [y]))
        posCoords = foldr (\next -> \prev -> cartProd next prev) [[]] ranges
    in map (\coords -> Position (coords!!0) (coords!!1) (coords!!2) 0.0) posCoords
         
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
            neighborIsUpper = isUpperSide side  
            f first = if neighborIsUpper && first then neighborVal else thisVal
        in if direction == directionFromCenter side 
            then
                let first = f True
                    second = f False
                in case (first, second) of
                    (Constant _, Constant c1) -> 
                        distributeMultiply [ first , Constant ( (-1)* c1) ] (1/interval)
                    _ -> distributeMultiply 
                            ( 
                                approximateDerivative vs 
                                    (if neighborIsUpper then neighborVal else thisVal) 
                                    (if neighborIsUpper then neighbor else position)  
                                ++ distributeMultiply
                                    (approximateDerivative vs 
                                        (if not neighborIsUpper then neighborVal else thisVal) 
                                        (if not neighborIsUpper then neighbor else position))
                                    (-1)
                            ) 
                            (1/interval)
            else 
                let proxyNeighbor p t = func (offsetPosition p $ t $ boundaryPair direction) Center 
                    thisProxyNeighborA = proxyNeighbor position fst
                    thisProxyNeighborB = proxyNeighbor position snd
                    neighborProxyNeighborA = proxyNeighbor neighbor fst
                    neighborProxyNeighborB = proxyNeighbor neighbor snd 
                in case ( thisProxyNeighborA, thisProxyNeighborB,neighborProxyNeighborA,neighborProxyNeighborB) of
                    (Constant thisA, Constant thisB, Constant nA, Constant nB) ->
                        distributeMultiply [ Constant (average [thisA,nA])
                            , Constant ((-1)* average [thisB,nB]) ] 
                            (1/ sl vs)    
                    _->undefined -- this is just too much work 
    _ -> []

solveUnknown::(Fractional a)=> ValSet a->Equation (Term a)->Position->a
solveUnknown vs (Equation l r unknownProperty) position= 
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
        U

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
    let 
    in concatMap (multTerm m) terms
        
multTerm m term  = case term of
    Constant c -> [Constant (c * m)]
    Unknown u -> [Unknown (u * m)]
    SubExpression s -> distributeMultiply ( getTerms s  ) m
    Derivative direc func side-> 
        let modf x s =  SubExpression $ Expression $ multTerm m $ func x s
        in [Derivative direc modf side]
                
integSurface:: (Num a)=> (Side->Term a) -> Position -> Direction -> Reader (ValSet a) [Term a]
integSurface f position direction = do
    vs <- ask
    return $ 
        let sides = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f 
                        else (\x ->
                            let [res] = multTerm (-1) x
                            in res  ).f 
                    term = modf s
                    sideAreaVal = runReader (sideArea s position) vs
                    nonDerivConstructor = case (direcDimenType direction,isUpper) of
                        (Temporal,True) -> Unknown
                        _ -> Constant 
                in case term of
                    Derivative _ subf _ -> 
                        SubExpression $ Expression $ distributeMultiply [subf position s] sideAreaVal
                    Constant c ->  nonDerivConstructor $ c * sideAreaVal
                    Unknown u ->  nonDerivConstructor $ u * sideAreaVal
        in [value (fst sides) True , value (snd sides) False]       
       
integSingleTerm::  Term Double -> DimensionType -> Position -> Reader (ValSet Double) [Term Double]
integSingleTerm term dimetype cellposition =  do
    vs <- ask
    return $
        let nonDerivativeAnswer = distributeMultiply [term] $ runReader (volumeOrInterval dimetype cellposition) vs 
        in case term of     
            Derivative direction func _ ->  
                if dimetype == direcDimenType direction 
                    then runReader (integSurface ( func cellposition ) cellposition direction) vs
                    else nonDerivativeAnswer
            _ -> nonDerivativeAnswer  

integ::  ValSet Double -> DimensionType -> [Term Double] ->Position -> [Term Double]
integ vs dimetype terms cellposition = case terms of
    [] -> []
    (x:xs) -> runReader (integSingleTerm x dimetype cellposition) vs ++ integ vs dimetype xs cellposition   

dd_:: (Fractional a) =>Property->Property-> Direction -> Reader (ValSet a) (Term a)
dd_ p1 p2 direction =do
    d<-ask 
    return $ Derivative direction (\x-> \s -> Constant $ runReader (prop p1 x s) d * runReader (prop p2 x s) d) Center

d_:: (Fractional a) => Property-> Direction -> Reader (ValSet a) (Term a)
d_ property direction =do
    d<-ask 
    return $ Derivative direction (\x-> \s -> Constant $ runReader (prop property x s) d) Center       
      
drho_dt = d_ Density Time
drhodu_dt =dd_ Density U Time   
drhodu_dx =dd_ Density U X
drhodv_dy = dd_ Density V Y
drhodw_dz = dd_ Density W Z   
 
continuity:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            [ integrate Temporal [runReader drho_dt env] >>= integrate Spatial , 
              integrate Spatial [runReader drhodu_dx env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodv_dy env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodw_dz env] >>= integrate Temporal ] 
            [\_ -> [Constant 0]] 
            Density
    
uMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
uMomentum = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            [ integrate Temporal [runReader drho_dt env] >>= integrate Spatial , 
              integrate Spatial [runReader drhodu_dx env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodv_dy env] >>= integrate Temporal, 
              integrate Spatial [runReader drhodw_dz env] >>= integrate Temporal ] 
            [\_ -> [Constant 0]] 
            U

vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
vMomentum = undefined    

wMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
wMomentum = undefined    

energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
energy = undefined    

gasLaw:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
gasLaw= undefined        
    
getDiscEqInstance:: Equation (Position -> [Term a]) -> Position -> Equation (Term a)
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) l) (concatMap (\t -> t pos) r) up
    
applyDiffEq :: (Fractional a)=>ValSet a -> Equation (Position -> [Term a]) -> ValSet a    
applyDiffEq (ValSet p v av sl) eq =
    let newVals = foldr
            (\pos -> \dict -> 
                let subDict = fromJust $ Map.lookup pos dict
                    discEquation=getDiscEqInstance eq pos 
                    solvedProperty = unknownProperty discEquation
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

testContEq:: Equation (Term Double)
testContEq = getDiscEqInstance ( runReader continuity initialGrid) testPosition 

ezcontinuity:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
ezcontinuity = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            [ integrate Temporal [runReader drho_dt env] >>= integrate Spatial,
              integrate Spatial [runReader drhodu_dx env] >>= integrate Temporal] 
            [\_ -> [Constant 0]]
            Density 
            
writeTermsOrig:: (Num a, Show a)=> [Term a] -> String
writeTermsOrig terms =
    let writeTerm t prev = prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) ++ " + "
            Derivative d f side -> "derivative( " ++(show d) ++ " " ++(show  side)++" )"
    in foldr writeTerm " " (terms |> reverse)  

writeTerms terms =
    let (x:y:xs) = writeTermsOrig terms |> reverse
    in xs |> reverse
  
testPosition =   Position 8 3 8 0.0
    
main:: IO()
main = 
    putStrLn " here we go... "
    -- >>= (\_-> print ( solveUnknown initialGrid testEquation $ Position 0 0 0 0.0)) 
    -- >>= (\_ -> putStrLn $ writeTerms $ distributeMultiply testTerms 2)
    -- >>= (\_ -> putStrLn $ show $ runReader (prop U testPosition Center ) initialGrid)
    >>= (\_ -> putStrLn $ writeTerms $ rhs testContEq)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs testContEq)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid testContEq testPosition)
    

