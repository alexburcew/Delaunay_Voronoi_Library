open Microsoft.Research.Science.Data.Imperative

type DataSet = Microsoft.Research.Science.Data.DataSet
type Vertex = Delaunay_Voronoi_Library.Vertex

type Hash_tbl<'Tkey,'Tval> = System.Collections.Generic.Dictionary<'Tkey,'Tval>

let mon_idx = 0
let data_var_name="temp"
let mv = -9999.0

let assert_ext cond message = System.Diagnostics.Debug.Assert(cond,message)


[<EntryPoint>]
let main argv =     
    use inputDs = DataSet.Open "GHCN_monthly.nc?openMode=readOnly"
    let N = inputDs.Dimensions.["station"].Length
    
    let lats = inputDs.GetData<single[]>("lat") |> Array.map (fun a -> float(a))
    let lons = inputDs.GetData<single[]>("lon") |> Array.map (fun a -> float(a))
    let ref_vals = inputDs.GetData<float[]>(sprintf "%s_mean" data_var_name,DataSet.ReduceDim(0),DataSet.Range(0,N-1))
    let ref_sigmas = inputDs.GetData<float[]>(sprintf "%s_sigma" data_var_name,DataSet.ReduceDim(0),DataSet.Range(0,N-1))

    let non_mv_indeces = [| for i in 0..N-1 do if ref_vals.[i]<>mv && ref_sigmas.[i]<>mv then yield i|]
    let non_mv_ind_set = Set<int>(non_mv_indeces)



    //excluding MVs    
    let excluder s =
        Seq.mapi (fun i e -> if non_mv_ind_set.Contains(i) then Some(e) else None) s |> Seq.choose (fun e -> e)

    let N = non_mv_indeces.Length
    let lats = excluder lats |> Seq.toArray
    let lons = excluder lons |> Seq.toArray
    let ref_vals = excluder ref_vals |> Seq.toArray
    let ref_sigmas = excluder ref_sigmas |> Seq.toArray

    let comp_vals = Array.zeroCreate<float> N
    let comp_sigmas = Array.zeroCreate<float> N

    let vertex_comarer = VertexComparer.Comparer()

    let variogram dist = 
        let s = 648.45
        let n = 45.46
        let r = 1492.2
        let alpha = 0.5
        let res = (s-n)*(1.0-exp(-dist*dist/(alpha*r*r)))+n
        res

    let covariogram dist =
        let s =32.215
        let l = 1059.6
        let res = s*s*exp(-abs(dist)/l)
        res

    let correlogram dist =
        let s = 0.90417
        let n = -0.0023
        let r = 2254.7
        let alpha = 0.5
        (s-n)*(exp(-dist*dist/(alpha*r*r)))+n

    let dist (v1:Vertex) (v2:Vertex)=
        let earthRadius = 6371.0
        let xdiff = earthRadius*(v1.X-v2.X)
        let ydiff = earthRadius*(v1.Y-v2.Y)
        let zdiff = earthRadius*(v1.Z-v2.Z)
        sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff)

    let hashTable = Hash_tbl<Vertex,float>(vertex_comarer)
    let vertices = [| for j in 0..N-1 do yield Vertex(lons.[j],lats.[j],ref_vals.[j]) |]
    for j in 0..N-1 do hashTable.Add(vertices.[j],ref_sigmas.[j])
    
    assert(hashTable.Count = N)

    let cov_0 = 22.0*22.0+covariogram(0.0)

    let computator i =
        printf "."
        if i%1000=0 then printf "%d" i
        let remove_idx idx data =
            Seq.append (Seq.take idx data) (Seq.skip (idx+1) data) |> Seq.toArray
        
        let ref_val = ref_vals.[i]
        let ref_sigma = ref_sigmas.[i]
        
        //cutting out ref values
        let eff_lats= remove_idx i lats
        let eff_lons = remove_idx i lons
        let eff_vals = remove_idx i ref_vals
        let eff_sigmas = remove_idx i ref_sigmas
        let vertices = [| for j in 0..N-2 do yield Vertex(eff_lons.[j],eff_lats.[j],eff_vals.[j]) |]
        let precomp = Delaunay_Voronoi_Library.Delaunay_Voronoi(System.Collections.Generic.List(vertices),false)            
        let res = precomp.NatNearestInterpolation(lons.[i],lats.[i],false,false)
        
        let lambdas = Array.map snd res.LambdasArray
        let verts = Array.map fst res.LambdasArray
        let m = verts.Length
                
        assert_ext (let b = Array.sum lambdas in abs(b-1.0)<=1e-10) "lambdas not normalized"
        assert_ext (Array.forall (fun a -> a>=0.0 && a<=1.0) lambdas) "lambda is out of 0.0 - 1.0"

        let sigma =
            let mutable acc = 0.0
            let mutable acc2 = 0.0
            for i in 0..m-1 do
                for j in 0..m-1 do
                    let cov = covariogram (dist verts.[i] verts.[j])
                    assert(cov>=0.0)
                    acc <- acc + lambdas.[i]*lambdas.[j]*cov
                        //if i<>j then
                        //else                            
                        //     hashTable.[verts.[i]]
                let cov2 = covariogram(dist verts.[i] res)
                assert(cov2>=0.0  && cov2<1040.0)
                acc2 <- acc2 + lambdas.[i]*cov2
            sqrt(cov_0 + acc + (-2.0*acc2))
                
        comp_vals.[i] <- res.Value;
        comp_sigmas.[i] <- sigma                    

    let tasks = [| for i in 0..N-1 do yield async { computator i; return() } |]

    Async.Parallel tasks |> Async.RunSynchronously
    

    use outputDs = DataSet.Open "out.csv?openMode=create&appendMetadata=false"
    outputDs.IsAutocommitEnabled <- false
    outputDs.AddVariable<float>("RefVal",ref_vals,"i") |> ignore
    outputDs.AddVariable<float>("RefSigma",ref_sigmas,"i") |> ignore
    outputDs.AddVariable<float>("Lat",lats,"i") |> ignore
    outputDs.AddVariable<float>("Lon",lons,"i") |> ignore
    outputDs.AddVariable<float>("CompVal",comp_vals,"i") |> ignore
    outputDs.AddVariable<float>("CompSigma",comp_sigmas,"i") |> ignore
    outputDs.Commit()

    0 // return an integer exit code
