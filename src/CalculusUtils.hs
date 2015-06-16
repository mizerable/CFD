module CalculusUtils where

import SolutionDomain
import Data.Maybe
import qualified Data.Map.Strict as Map
import Control.Monad.Reader as Reader
import Data.List


timeStep :: (Num a,Fractional a) => a            
timeStep = 0.01

specificHeatCv :: (Num a) => a
specificHeatCv = 15

-- sideArea:: (Num a, Fractional a)=>Side -> Position -> a
sideArea s (Position x y z _) = case s of 
    Now -> 1
    Prev -> 1
    _ -> fromJust $! Map.lookup (Position x y z 0) (areaVal $! initialGrid)  >>= Map.lookup s    

-- sideLength:: (Num a, Fractional a) => Direction -> Position ->  a
sideLength d (Position x y z _) = case d of 
    Time -> timeStep
    _ -> fromJust $! Map.lookup (Position x y z 0) (sideLen $! initialGrid) >>= Map.lookup d   

boundaryPair:: Direction -> (Side,Side)
boundaryPair d = case d of 
     X -> (East,West)
     Y -> (North,South)
     Z -> (Top,Bottom)
     Time -> (Now,Now) -- we don't use the previous time step, actually it's more like (future ,now)  but we don't have a future
         
direcDimenType:: Direction -> DimensionType
direcDimenType direc = case direc of
    Time -> Temporal
    _ -> Spatial

--volumeOrInterval:: (Num a, Fractional a ) => DimensionType -> Position ->a  
volumeOrInterval dimetype position = case dimetype of
    Temporal -> timeStep
    Spatial -> foldl' (\p d  -> p * sideLength d position ) 1 $! enumFrom X 

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

-- approximateDerivative::(Num a, Fractional a)=>  Term a -> Position-> Term a
approximateDerivative deriv position= case deriv of 
    (Derivative direction func side m) ->
        let neighbor = offsetPosition position side
            sl =  sideLength direction position
            sln = sideLength direction neighbor
            interval = average [sl, sln ]
            thisVal = func position Center 
            neighborVal = func neighbor Center
            neighborIsUpper = isUpperSide side  
            f first = case (first, neighborIsUpper) of
                (True, True) -> neighborVal
                (False, False) -> neighborVal
                _-> thisVal
            mult = m position Center
        in case (f True, f False) of
            (Constant c1 , Constant c2) -> Constant $ (c1-c2)*mult/interval 
            _ -> error "can't approx >1 order derivs. deriv must produce constants" 
    _ -> error "can't approx something that isn't a deriv"

--solveUnknown::(Fractional a)=> Equation (Term a)->Position->a
solveUnknown (Equation l r _) position= 
    let sumUnknown p n =  p + case n of
            Unknown u-> u
            SubExpression s -> sumExpression sumUnknown $! getTerms s
            _ -> 0
        sumConstants p n =  p + case n of
            Constant c-> c
            Derivative {}-> sumExpression sumConstants $!
                [approximateDerivative n $! position]  
            SubExpression s -> sumExpression sumConstants $! getTerms s
            _ -> 0
        sumExpression s = foldl' s 0
        lhsUnknown = sumExpression sumUnknown l
        rhsUnknown = sumExpression sumUnknown r
        lhsConstants = sumExpression sumConstants l
        rhsConstants = sumExpression sumConstants r
    in (rhsConstants - lhsConstants)/(lhsUnknown-rhsUnknown)

distributeMultiply::(Num a)=> [Term a]->a->[Term a]
distributeMultiply terms m = concatMap (multTerm m) terms
        
multTerm:: (Num a)=> a -> Term a -> [Term a]
multTerm m term  = case term of
    SubExpression s -> distributeMultiply ( getTerms s  ) m
    Derivative direc func side multF->  
        let modf x s = fmap (*m) (func x s)   
        in [Derivative direc modf side multF]
    _-> [fmap (*m) term]
                
--integSurface:: (Num a,Fractional a, RealFloat a)=> (Side->Term a) -> Position -> Direction ->Property-> Reader (ValSet a) [Term a]
integSurface f position direction unknownProp= do
    vs <- ask
    return $ 
        let sides = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f 
                        else (\x -> let [res] = multTerm (-1) x
                                in res  ).f 
                    term = modf s
                    sideAreaVal = sideArea s position
                    isUnknown = (direcDimenType direction,isUpper) == (Temporal,True) 
                in case term of
                    Derivative d subf _ m-> head $ multTerm sideAreaVal $ Derivative d subf s m
                    SubExpression (Expression expr) -> SubExpression $ Expression $ distributeMultiply expr sideAreaVal
                    _-> if isUnknown 
                        then 
                            let constantVal = fmap
                                    (* (sideAreaVal 
                                        / prop unknownProp position s vs)) 
                                    term
                            in Unknown $ if isNaN (val constantVal) 
                                then 1 else val constantVal      
                        else fmap (* sideAreaVal) term
        in [value (fst sides) True , value (snd sides) False]       
       
--integSingleTerm::  (Num a, Fractional a, RealFloat a ) =>  Term a -> DimensionType -> Position -> Property ->Reader (ValSet a) [Term a]
integSingleTerm term dimetype cellposition unknownProp=  do
    vs <- ask
    return $
        let nonDerivAnswer = case term of 
                SubExpression _ -> error "can't integrate a subexpression as a single term"
                _ -> multTerm (volumeOrInterval dimetype cellposition) term
                --Unknown c -> multTerm (runReader (volumeOrInterval dimetype cellposition) vs) term -- error( "Unknown " ++ show c)
                --Constant c -> multTerm (runReader (volumeOrInterval dimetype cellposition) vs) term --   error ("Constant " ++ show c)
                --(Derivative d f c m) -> 
                    --multTerm (runReader (volumeOrInterval dimetype cellposition) vs) (Derivative X (\_ _ -> Constant 1) Center (\_ _ -> 1))
                    --multTerm (runReader (volumeOrInterval dimetype cellposition) vs) term
                    --error ("Derivative " ++ (show d) ++ " " ++ (show c ) )   
        in case term of     
            Derivative direction func _ _->
                if direcDimenType direction == dimetype  
                then runReader 
                        (integSurface ( func cellposition ) cellposition direction unknownProp)
                        vs
                else nonDerivAnswer
            _ -> nonDerivAnswer    

--integ::  (Num a, Fractional a, RealFloat a ) => ValSet a-> DimensionType -> [Term a] ->Position ->Reader Property [Term a]
integ vs dimetype terms cellposition = do
    unknownProp <- ask
    return $ case terms of
        [] -> []
        (x:xs) -> runReader 
            (integSingleTerm x dimetype cellposition unknownProp) vs 
            ++ runReader (integ vs dimetype xs cellposition) unknownProp   

--d_:: (Fractional a) => [Property]-> Direction -> Reader (ValSet a) (Term a)
d_ properties = df_ properties 1 

-- df_ :: (Num a,Fractional a) => [Property]-> a ->Direction -> Reader (ValSet a) (Term a)
df_ :: [Property]-> Double ->Direction -> Reader (ValSet Double) (Term Double)      
df_ properties factor = dfm_ properties factor (\_ _ -> 1)         

dfm_ properties factor multF direction = do
    d<-ask 
    return $ Derivative direction 
        (\x s -> Constant $ runReader ( multProps properties x s ) d   * factor )
        Center
        multF       
      
dd_ :: (Num a,Fractional a)=> Reader (ValSet a) (Term a) -> Direction -> Reader (ValSet a) (Term a)
dd_ inner  = ddf_ inner 1 

ddf_ ::(Num a,Fractional a)=> Reader (ValSet a) (Term a) -> a-> Direction -> Reader (ValSet a) (Term a)
ddf_ inner factor = ddfm_ inner factor (\_ _ -> return 1)

ddfm_ inner factor multF direction = do
    d<- ask
    return $ Derivative 
        direction 
        (\_ _ -> head $! multTerm factor $! runReader inner d  ) 
        Center 
        (\pos side -> runReader (multF pos side) d)
       
--drho_dt:: (Fractional a) => Reader (ValSet a) (Term a)        
drho_dt = d_ [Density] Time

--dp_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dx = d_ [Pressure] X

--dp_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dy = d_ [Pressure] Y

--dp_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dp_dz = d_ [Pressure] Z

--drhou_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhou_dt =d_ [Density, U] Time   

--drhov_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhov_dt =d_ [Density, V] Time

--drhow_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhow_dt =d_ [Density, W] Time

--drhoT_dt:: (Fractional a) => Reader (ValSet a) (Term a)
drhoT_dt = df_ [Density, Temperature] specificHeatCv Time 

--dmewu_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dx = d_ [Mew, U] X

--dmewu_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dy = d_ [Mew, U] Y

--dmewu_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmewu_dz = d_ [Mew, U] Z

--dmewv_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dx = d_ [Mew, V] X

--dmeww_dx:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dx = d_ [Mew, W] X

--dmewv_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dy = d_ [Mew,V] Y

--dmewv_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmewv_dz = d_ [Mew,V] Z

--dmeww_dy:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dy = d_ [Mew,W] Y

--dmeww_dz:: (Fractional a) => Reader (ValSet a) (Term a)
dmeww_dz = d_ [Mew,W] Z

--divergence:: (Num a, Fractional a) => [Reader (ValSet a) (Term a)] -> [Reader (ValSet a) (Term a)]
divergence vector = map ( uncurry dd_) $! zip vector (enumFrom X)  
    
-- gradient:: (Num a, Fractional a) => [Property] -> a-> [Reader (ValSet a) (Term a)]
gradient properties constantFactor = map (df_ properties constantFactor) $! enumFrom X  

integrateTerms integrate env =map (\term -> integrateOrder integrate Spatial Temporal [runReader term env] )  

divergenceWithProps properties = divergenceWithPropsFactor properties 1  

divergenceWithPropsFactor properties constantFactor = 
    map (\x -> df_ (properties++[d2C x]) constantFactor x ) (enumFrom X)
 
divGrad properties constantFactor = divergence $ gradient properties constantFactor  

integrateOrder integrate i1 i2 term = integrate i1 term >>= integrate i2
        
-- multProps ::(Num a, Fractional a)=> [Property] ->Position ->Side -> Reader (ValSet a) a
multProps = 
    foldl' 
        (\prev next pos side->
            do
                vs <- ask
                return $ runReader (prev pos side) vs *  prop next pos side vs)         
        (\_ _ -> return 1) 
                
-- squareDerivative:: (Num a, Fractional a) => [Property] -> a->Direction -> [Reader (ValSet a) (Term a)]        
squareDerivative properties constantFactor direction = 
    let foldedProps = multProps properties
    in [ ddf_ (d_ (properties++properties) direction) (constantFactor/2) direction
        ,ddfm_ (d_ properties direction) (constantFactor * (-1)) foldedProps direction] 

--pairedMultipliedDerivatives :: (Num a, Fractional a) =>  [Property] -> [Property] -> Direction -> Direction -> [Reader (ValSet a) (Term a)]
pairedMultipliedDerivatives props1 props2 dir1 dir2 = 
    let p1 = multProps props1 
        p2 = multProps props2
    in [dd_ (d_ (props1++props2) dir1) dir2,
        ddfm_ (d_ props1 dir1) (-1) p2 dir2,
        ddfm_ (d_ props2 dir1) (-1) p1 dir2]

integUnknown env unknownProp dimetype terms cellposition = 
    runReader (integ env dimetype terms cellposition)unknownProp



