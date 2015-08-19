{-# LANGUAGE ScopedTypeVariables, RankNTypes , KindSignatures , FlexibleContexts #-}
module CalculusUtils where

import SolutionDomain
import Control.Monad.Reader as Reader
import Data.List
         
direcDimenType:: Direction -> DimensionType
direcDimenType direc = case direc of
    Time -> Temporal
    _ -> Spatial

--volumeOrInterval:: (Num a, Fractional a ) => DimensionType -> Position ->a  
volumeOrInterval dimetype position vs = case dimetype of
    Temporal -> timeStep
    Spatial -> foldl' (\p d  -> p * sideLength d position vs) 1 $! enumFrom X 

c2D:: Property -> Direction
c2D c = case c of
    U -> X
    V -> Y
    W -> Z
    _ -> error "other properties don't have an associated direction"
    
d2C:: Direction -> Property
d2C d = case d of 
    X -> U
    Y -> V
    Z -> W
    _ -> error "other directions don't have an associated property"

--solveUnknown::(Fractional a)=> Equation (SchemeType ->Term a)->Position->a
solveUnknown (Equation l r _) position vs= 
    let sumUnknown p n =  p + case n of
            Unknown u-> u
            SubExpression s -> sumExpression sumUnknown $! getTerms s
            _ -> 0
        sumConstants p n =  p + case n of
            Constant c-> c
            Derivative {}-> sumExpression sumConstants [approximateDerivative n position vs]  
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
        
--multTerm:: (Num a)=> a -> Term a -> [Term a]
multTerm m term  = case term of
    SubExpression s -> distributeMultiply ( getTerms s  ) m
    Derivative direc func side multF->  
        let modf st x s = fmap (*m) (func st x s)   
        in [Derivative direc modf side multF]
    _-> [fmap (*m) term]
                
--integSurface:: (Num a,Fractional a, RealFloat a)=> (Side->Term a) -> Position -> Direction ->Property-> Reader (ValSet a) [Term a]
integSurface f position direction unknownProp= do
    vs <- ask
    return $ 
        let (s1, s2) = boundaryPair direction 
            value s isUpper =
                let modf = if isUpper then f 
                        else (\x -> let [res] = multTerm (-1) x
                                in res  ).f 
                    term = modf s
                    sideAreaVal = sideArea s position vs
                    isUnknown = (direcDimenType direction,isUpper) == (Temporal,True) 
                in case term of
                    Derivative d subf _ m-> head $ multTerm sideAreaVal $ Derivative d subf s m
                    SubExpression (Expression expr) -> SubExpression $ Expression $ distributeMultiply expr sideAreaVal
                    _-> if isUnknown 
                        then 
                            let constantVal = fmap
                                    (* (sideAreaVal 
                                        / prop Directional unknownProp position s vs)) 
                                    term
                            in Unknown $ if isNaN (val constantVal) 
                                then 1 else val constantVal      
                        else fmap (* sideAreaVal) term
        in [value s1 True , value s2 False]       
       
integSingleTerm :: forall (m :: * -> *).
                         MonadReader (ValSet Double) m =>
                         Term Double
                         -> DimensionType -> Position -> Property -> m [Term Double]
integSingleTerm term dimetype cellposition unknownProp=  do
    vs <- ask
    return $
        let nonDerivAnswer = case term of 
                SubExpression _ -> error "can't integrate a subexpression as a single term"
                _ -> multTerm (volumeOrInterval dimetype cellposition vs) term   
        in case term of     
            Derivative direction func _ _->
                if direcDimenType direction == dimetype  
                then runReader 
                        (integSurface ( func Directional cellposition ) cellposition direction unknownProp)
                        vs
                else nonDerivAnswer
            _ -> nonDerivAnswer    

integ vs dimetype terms cellposition = do
    unknownProp <- ask
    return $ case terms of
        [] -> []
        (x:xs) -> runReader 
            (integSingleTerm x dimetype cellposition unknownProp) vs 
            ++ runReader (integ vs dimetype xs cellposition) unknownProp   

--d_::[Property]-> Direction -> Reader (ValSet Double) (Term Double)
d_ properties = df_ properties 1 

-- df_ :: (Num a,Fractional a) => [Property]-> a ->Direction -> Reader (ValSet a) (Term a)
df_ :: [Property]-> Double ->Direction -> Reader (ValSet Double) (Term Double)      
df_ properties factor = dfm_ properties factor (\_ _ _ -> 1)         

dfm_ :: forall (m :: * -> *).
          MonadReader (ValSet Double) m =>
          [Property]
          -> Double
          -> (SchemeType -> Position -> Side -> Double)
          -> Direction
          -> m (Term Double)
dfm_ properties factor multF direction = do
    d<-ask 
    return $ Derivative direction 
        (\st x s -> Constant $ runReader (  multProps properties st x s ) d   * factor )
        Center
        multF       
      
--dd_ :: (Num a,Fractional a)=> Reader (ValSet a) (Term a) -> Direction -> Reader (ValSet a) (Term a)
dd_ inner  = ddf_ inner 1 

--ddf_ ::(Num a,Fractional a)=> Reader (ValSet a) (Term a) -> a-> Direction -> Reader (ValSet a) (Term a)
ddf_ inner factor = ddfm_ inner factor (\_ _ _ -> return 1)

ddfm_ inner factor multF direction = do
    d<- ask
    return $ Derivative 
        direction 
        (\_ _ _-> head $! multTerm factor $! runReader inner d  ) 
        Center 
        (\st pos side -> runReader (multF st pos side) d)
       
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

drhodye_dt =d_ [Density, Dye] Time

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

divergence :: forall (m :: * -> *) a r.
                    (MonadReader r m, Num a) =>
                    [Reader r (Term a)] -> [m (Term a)]
divergence vector = map ( uncurry dd_) $! zip vector (enumFrom X)  
    
gradient :: [Property] -> Double -> [Reader (ValSet Double) (Term Double)]
gradient properties constantFactor = map (df_ properties constantFactor) $! enumFrom X  

integrateTerms :: forall (m :: * -> *) t r.
                    Monad m =>
                    (DimensionType -> [t] -> m [t]) -> r -> [Reader r t] -> [m [t]]
integrateTerms integrate env =map (\term -> integrateOrder integrate Spatial Temporal [runReader term env] )  

divergenceWithProps :: [Property]-> [Reader (ValSet Double) (Term Double)]
divergenceWithProps properties = divergenceWithPropsFactor properties 1  

divergenceWithPropsFactor :: [Property] -> Double -> [Reader (ValSet Double) (Term Double)]
divergenceWithPropsFactor properties constantFactor = 
    map (\x -> df_ (properties++[d2C x]) constantFactor x ) (enumFrom X)
 
divGrad :: forall (m :: * -> *).
                 MonadReader (ValSet Double) m =>
                 [Property] -> Double -> [m (Term Double)]
divGrad properties constantFactor = divergence $ gradient properties constantFactor  

integrateOrder :: forall (m :: * -> *) b t.
                        Monad m =>
                        (t -> b -> m b) -> t -> t -> b -> m b
integrateOrder integrate i1 i2 term = integrate i1 term >>= integrate i2
        
multProps props schemeType= 
    let func = prop schemeType
    in foldl' 
        (\prev next pos side->
            do
                vs <- ask
                return $ runReader (prev pos side) vs *  func next pos side vs)          
        (\_ _ -> return 1)
        props  
      
squareDerivative :: forall (m :: * -> *).
                          MonadReader (ValSet Double) m =>
                          [Property] -> Double -> Direction -> [m (Term Double)]
squareDerivative properties constantFactor direction = 
    let foldedProps = multProps properties
    in [ ddf_ (d_ (properties++properties) direction) (constantFactor/2) direction
        ,ddfm_ (d_ properties direction) (constantFactor * (-1)) foldedProps direction] 

pairedMultipliedDerivatives :: forall (m :: * -> *).
                                     MonadReader (ValSet Double) m =>
                                     [Property]
                                     -> [Property] -> Direction -> Direction -> [m (Term Double)]
pairedMultipliedDerivatives props1 props2 dir1 dir2 = 
    let p1 = multProps props1 
        p2 = multProps props2
    in [dd_ (d_ (props1++props2) dir1) dir2,
        ddfm_ (d_ props1 dir1) (-1) p2 dir2,
        ddfm_ (d_ props2 dir1) (-1) p1 dir2]

integUnknown :: ValSet Double
                  -> Property
                  -> DimensionType
                  -> [Term Double]
                  -> Position
                  -> [Term Double]
integUnknown env unknownProp dimetype terms cellposition = 
    runReader (integ env dimetype terms cellposition)unknownProp



