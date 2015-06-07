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
data Property = U | V | W | Density | Temperature | Mew | Pressure deriving (Ord,Eq,Enum,Show)
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
    | Derivative { denom::Direction ,function:: Position->Side->Term a, centered::Side, 
        multiplier:: Position->Side-> a }

instance Functor Term  where
    fmap f x = case x of
         Constant c -> Constant $ f c
         Unknown u -> Unknown $ f u
         _ -> undefined    
            
timeStep::Double
timeStep = 0.0001 

specificHeatCv :: Double
specificHeatCv = 15

heatConductivityK:: Double
heatConductivityK = 0.1

maxPos:: Direction -> Int
maxPos d = case d of 
    X -> 12
    Y -> 12
    Z -> 12
    Time -> undefined

getPositionComponent:: Position -> Direction -> Int
getPositionComponent (Position x y z _) d = case d of 
    X -> x
    Y -> y
    Z -> z
    Time -> undefined

c2D:: Property -> Direction
c2D c = case c of
    U -> X
    V -> Y
    W -> Z
    _ -> undefined
    
d2C:: Direction -> Property
d2C d = case d of 
    X -> U
    Y -> V
    Z -> W
    _ -> undefined

modifyPositionComponent::Position -> Direction -> Int -> Position
modifyPositionComponent (Position x y z t) direction amt= case direction of 
    X -> Position amt y z t
    Y -> Position x amt z t 
    Z -> Position x y amt t
    Time -> undefined -- idk what to do here, time is a double and it won't let me add
    
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
            $ maxOrMin boundary $ getPositionComponent position direction + offsetAmount   
        
prop::(Num a, Fractional a)=> Property->Position->Side-> Reader (ValSet a) a
prop property position side = do
    (ValSet _ v _ _) <- ask 
    return $ 
        let neighbor = offsetPosition position side
            getVal p = fromJust $ Map.lookup p v >>= Map.lookup property  
        in average [getVal position,getVal neighbor]

initialGrid::(Num a,Fractional a) => ValSet a
initialGrid= 
    let p = makePositions
        vMap = foldr (\next -> \prev -> Map.insert next 
            (case next of 
                U-> 1.153
                V-> 0
                W-> 0
                _-> 1
                ) 
            prev) Map.empty (enumFrom U)
        avMap = foldr (\next -> \prev -> Map.insert next 1.1 prev) Map.empty (enumFrom East)
        slMap = foldr (\next -> \prev -> Map.insert next 1.12 prev) Map.empty (enumFrom X)
        v = foldr (\next -> \prev -> Map.insert next vMap prev) Map.empty makePositions
        av = foldr (\next -> \prev -> Map.insert next avMap prev) Map.empty makePositions
        sl = foldr (\next -> \prev -> Map.insert next slMap prev) Map.empty makePositions
    in ValSet p v av sl
    
cartProd:: [[a]] -> [[a]] -> [[a]]
cartProd xs ys = [ x ++ y | x <- xs, y <- ys]

makePositions::[Position]
makePositions = 
    let ranges = enumFrom X |> map maxPos |> map (\x -> [1..x] |> map (\y -> [y]))
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
    Center -> undefined 

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
    (Derivative direction func side multiplier) ->
        let neighbor = offsetPosition position side
            sl =  runReader (sideLength direction position)
            sln = runReader (sideLength direction neighbor)
            interval = average [sl vs, sln vs]
            thisVal = func position Center 
            neighborVal = func neighbor Center
            neighborIsUpper = isUpperSide side  
            f first = if neighborIsUpper && first then neighborVal else thisVal
            mult = multiplier position Center
        in if direction == directionFromCenter side 
            then
                let first = f True
                    second = f False
                in case (first, second) of
                    (Constant _, Constant c1) -> 
                        distributeMultiply [ first , Constant ( (-1)* c1) ] (mult /interval)
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
                            (mult /interval)
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
                            (mult / sl vs)    
                    _->undefined -- this is just too much work 
    _ -> []

solveUnknown::(Fractional a)=> ValSet a->Equation (Term a)->Position->a
solveUnknown vs (Equation l r _) position= 
    let sumUnknown n p =  p + case n of
            Unknown u-> u
            SubExpression s -> sumExpression sumUnknown $ getTerms s
            _ -> 0
        sumConstants n p =  p + case n of
            Constant c-> c
            Derivative _ _ _ _-> sumExpression sumConstants $ approximateDerivative vs n position  
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
distributeMultiply terms m = concatMap (multTerm m) terms
        
multTerm:: (Num a)=> a -> Term a -> [Term a]
multTerm m term  = case term of
    SubExpression s -> distributeMultiply ( getTerms s  ) m
    Derivative direc func side multF-> 
        let modf x s =  SubExpression $ Expression $ multTerm m $ func x s
        in [Derivative direc modf side multF]
    _-> [fmap (\x-> x*m) term]
                
integSurface:: (Num a)=> (Side->Term a) -> Position -> Direction -> Reader (ValSet a) [Term a]
integSurface f position direction = do
    vs <- ask
    return $ 
        let sides = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f 
                        else (\x -> let [res] = multTerm (-1) x
                                in res  ).f 
                    term = modf s
                    sideAreaVal = runReader (sideArea s position) vs
                    isUnknown = (direcDimenType direction,isUpper) == (Temporal,True) 
                in case term of
                    Derivative _ subf _ _-> 
                        SubExpression $ Expression $ distributeMultiply [subf position s] sideAreaVal
                    SubExpression (Expression expr) -> SubExpression $ Expression $ distributeMultiply expr sideAreaVal
                    _-> if isUnknown 
                        then Unknown $ sideAreaVal 
                        else fmap (\x-> x * sideAreaVal) term
        in [value (fst sides) True , value (snd sides) False]       
       
integSingleTerm::  Term Double -> DimensionType -> Position -> Reader (ValSet Double) [Term Double]
integSingleTerm term dimetype cellposition =  do
    vs <- ask
    return $
        let nonDerivativeAnswer = distributeMultiply [term] $ runReader (volumeOrInterval dimetype cellposition) vs 
        in case term of     
            Derivative direction func _ _->  
                if dimetype == direcDimenType direction 
                    then runReader (integSurface ( func cellposition ) cellposition direction) vs
                    else nonDerivativeAnswer
            _ -> nonDerivativeAnswer  

integ::  ValSet Double -> DimensionType -> [Term Double] ->Position -> [Term Double]
integ vs dimetype terms cellposition = case terms of
    [] -> []
    (x:xs) -> runReader (integSingleTerm x dimetype cellposition) vs ++ integ vs dimetype xs cellposition   

d_:: (Fractional a) => [Property]-> Direction -> Reader (ValSet a) (Term a)
d_ properties direction = df_ properties 1 direction

df_ :: (Num a,Fractional a) => [Property]-> a ->Direction -> Reader (ValSet a) (Term a)      
df_ properties factor direction = dfm_ properties factor (\_ -> \_ -> 1) direction         

dfm_ properties factor multF direction = do
    d<-ask 
    return $ Derivative direction 
        (\x-> \s -> 
            foldr (\next -> \prev -> fmap ((*) ( runReader (prop next x s) d )) prev  |> 
                multTerm factor |> head) (Constant 1) properties )
        Center
        multF       
      
dd_ :: (Num a,Fractional a)=> Reader (ValSet a) (Term a) -> Direction -> Reader (ValSet a) (Term a)
dd_ inner direction = ddf_ inner 1 direction

ddf_ ::(Num a,Fractional a)=> Reader (ValSet a) (Term a) -> a-> Direction -> Reader (ValSet a) (Term a)
ddf_ inner factor direction = ddfm_ inner factor (\_ -> \_ -> 1) direction

ddfm_ inner factor multF direction = do
    d<- ask
    return $ Derivative direction (\x-> \s -> runReader inner d |> multTerm factor |> head ) Center multF
        
drho_dt = d_ [Density] Time

dp_dx = d_ [Pressure] X
dp_dy = d_ [Pressure] Y
dp_dz = d_ [Pressure] Z

drhou_dt =d_ [Density, U] Time   
drhov_dt =d_ [Density, V] Time
drhow_dt =d_ [Density, W] Time

drhoT_dt = df_ [Density, Temperature] specificHeatCv Time 

dmewu_dx = d_ [Mew, U] X
dmewu_dy = d_ [Mew, U] Y
dmewu_dz = d_ [Mew, U] Z

dmewv_dx = d_ [Mew, V] X
dmeww_dx = d_ [Mew, W] X

dmewv_dy = d_ [Mew,V] Y
dmewv_dz = d_ [Mew,V] Z
dmeww_dy = d_ [Mew,W] Y
dmeww_dz = d_ [Mew,W] Z

drhoTu_dx = df_ [Density,Temperature,U] specificHeatCv X
drhoTv_dy = df_ [Density,Temperature,V] specificHeatCv Y
drhoTw_dz = df_ [Density,Temperature,W] specificHeatCv Z

dpu_dx = d_ [Pressure,U] X
dpv_dy = d_ [Pressure,V] Y
dpw_dz = d_ [Pressure,W] Z

divergence:: (Num a, Fractional a) => [Reader (ValSet a) (Term a)] -> [Reader (ValSet a) (Term a)]
divergence vector = 
    let zipped = zip vector (enumFrom X) 
    in zipped |> map (\pair -> dd_ (fst pair) (snd pair) )
    
gradient:: (Num a, Fractional a) => [Property] -> a-> [Reader (ValSet a) (Term a)]
gradient properties constantFactor = enumFrom X |> map (\d-> df_ properties constantFactor d)

integrateTerms integrate env =map (\term -> integrateOrder integrate Spatial Temporal [runReader term env] )  

divergenceWithProps properties = divergenceWithPropsFactor properties 1  

divergenceWithPropsFactor properties constantFactor = 
    map (\x -> df_ (properties++[d2C x]) constantFactor x ) (enumFrom X)
 
divGrad properties constantFactor = (divergence $ gradient properties constantFactor)  

integrateOrder integrate i1 i2 term = integrate i1 term >>= integrate i2
        
squareDerivative properties direction = do
    vs <- ask
    return $ 
        let foldedProps = foldr 
                (\next -> \prev -> (\pos -> \side -> prev pos side * (runReader (prop next pos side) vs))) 
                (\_ -> \_ -> 1) 
                properties
        in [ddf_ (d_ (properties++properties) direction) (1/2) direction 
            ,ddfm_ (d_ properties direction) 1 foldedProps direction] 

continuity:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
continuity = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            ([ integrate Temporal [runReader drho_dt env] >>= integrate Spatial]
                ++ (integrateTerms integrate env $ divergenceWithProps [Density])  
            )
            [\_ -> [Constant 0]] 
            Density
    
uMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
uMomentum = do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            ([ integrate Temporal [runReader drhou_dt env] >>= integrate Spatial ]
                ++ (integrateTerms integrate env $ divergenceWithProps [Density, U])
                ++ [integrate Spatial [runReader dp_dx env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,U] 1 
                  ,divergence [ dmewu_dx, dmewv_dx, dmeww_dx ] 
                  ,map (\d-> ddf_ d (-2/3) X) (divergenceWithProps [Mew]) 
                ]
            )
            U

vMomentum::Reader (ValSet Double) (Equation (Position-> [Term Double]))
vMomentum = do
    env <- ask
    return $ let integrate = integ env
        in Equation
            ([ integrate Temporal [runReader drhov_dt env] >>= integrate Spatial ]
                ++ (integrateTerms integrate env $ divergenceWithProps [Density, V])
                ++ [integrate Spatial [runReader dp_dy env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,V] 1 
                  ,divergence [ dmewu_dy, dmewv_dy, dmeww_dy ] 
                  ,map (\d-> ddf_ d (-2/3) Y) (divergenceWithProps [Mew]) 
                ]
            )
            V   

wMomentum:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
wMomentum =  do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            ([ integrate Temporal [runReader drhow_dt env] >>= integrate Spatial ]
                ++ (integrateTerms integrate env $ divergenceWithProps [Density , W] ) 
                ++[integrate Spatial [runReader dp_dz env] >>= integrate Temporal ]
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Mew,W] 1 
                  ,divergence [ dmewu_dz, dmewv_dz, dmeww_dz ] 
                  ,map (\d-> ddf_ d (-2/3) Z) (divergenceWithProps [Mew]) 
                ]
            )
            W      

energy::Reader (ValSet Double) (Equation (Position-> [Term Double]))
energy =  do
    env <- ask
    return $ let integrate = integ env 
        in Equation
            ([ integrate Temporal [runReader drhoT_dt env] >>= integrate Spatial ]
                ++ (integrateTerms integrate env $ divergenceWithPropsFactor [Density , W] specificHeatCv ) 
                ++ (integrateTerms integrate env $ divergenceWithProps [Pressure] )
            ) 
            ( concatMap (integrateTerms integrate env) 
                [ divGrad [Temperature] heatConductivityK 
                   
                ]
            )
            Temperature   

gasLaw:: Reader (ValSet Double) (Equation (Position-> [Term Double]))
gasLaw= undefined        
    
getDiscEqInstance:: Equation (Position -> [Term a]) -> Position -> Equation (Term a)
getDiscEqInstance (Equation l r up) pos = Equation (concatMap (\t -> t pos) l) (concatMap (\t -> t pos) r) up
    
advancePositionTime :: Position -> Position
advancePositionTime (Position x y z t ) = Position x y z (t+timeStep)    
    
applyDiffEq :: (Fractional a)=>ValSet a -> Equation (Position -> [Term a]) -> ValSet a    
applyDiffEq (ValSet p v av sl) eq =
    let newVals = foldr
            (\pos -> \dict -> 
                let subDict = fromJust $ Map.lookup pos dict
                    discEquation=getDiscEqInstance eq pos 
                    solvedProperty = unknownProperty discEquation
                    newValue = solveUnknown (ValSet p v av sl) discEquation pos
                    newPos = advancePositionTime pos  
                in Map.insert newPos (Map.insert solvedProperty newValue subDict)dict )  
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

testEq:: Reader (ValSet Double ) (Equation (Position ->[Term Double])) -> Equation (Term Double)
testEq eq = getDiscEqInstance ( runReader eq initialGrid) testPosition 
            
writeTermsOrig:: (Num a, Show a)=> [Term a] -> String
writeTermsOrig terms =
    let writeTerm t prev = prev ++ case t of
            Unknown u -> show u ++ "X + "
            Constant c -> show c ++ " + "
            SubExpression s -> writeTerms (getTerms s) ++ " + "
            Derivative d _ side _-> "derivative( " ++ show d ++ " " ++ show  side ++" )"
    in foldr writeTerm " " (terms |> reverse)  

writeTerms terms =
    let (_:_:xs) = writeTermsOrig terms |> reverse
    in xs |> reverse
  
testPosition =   Position 8 3 9 0.0
    

    
main:: IO()
main = 
    putStrLn "starting ..... "
    >>= (\_-> print ( solveUnknown initialGrid testEquation $ Position 0 0 0 0.0)) 
    >>= (\_ -> putStrLn $ writeTerms $ distributeMultiply testTerms 2)
    >>= (\_ -> putStrLn $ show $ runReader (prop U testPosition Center ) initialGrid)
    >>= (\_ -> putStrLn " continuity ------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq continuity)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq continuity)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid (testEq continuity) testPosition)
    >>= (\_ -> putStrLn " U Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq uMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq uMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid (testEq uMomentum) testPosition)
    >>= (\_ -> putStrLn " V Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq vMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq vMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid (testEq vMomentum) testPosition)
    >>= (\_ -> putStrLn " W Momentum------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq wMomentum)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq wMomentum)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid (testEq wMomentum) testPosition)
    >>= (\_ -> putStrLn " ENERGY ------------ ")
    >>= (\_ -> putStrLn $ writeTerms $ rhs $ testEq energy)
    >>= (\_ -> putStrLn " = ")
    >>= (\_ -> putStrLn $ writeTerms $ lhs $ testEq energy)
    >>= (\_ -> putStrLn " solving... ")
    >>= (\_ -> putStrLn $ show $ solveUnknown initialGrid (testEq energy) testPosition)

