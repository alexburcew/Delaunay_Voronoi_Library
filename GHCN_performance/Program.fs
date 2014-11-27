open Microsoft.Research.Science.Data
open System.Diagnostics

type table = {
    lat: float[]
    lon: float[]
    value: float[]
    }

let trace = Trace.WriteLine

let subruns = 31

let predictGrid grid obs =
    let zipped = Array.zip3 obs.lat obs.lon obs.value
    let vertices = System.Collections.Generic.List(Seq.mapi (fun i x -> let lat,lon,v= x in Delaunay_Voronoi_Library.Vertex(lon,lat,v,0.0,i)) zipped)
    let sw = Stopwatch.StartNew()
    let delauney = new Delaunay_Voronoi_Library.Delaunay_Voronoi(vertices,false)
    sprintf "%A: delauley generated" sw.Elapsed |> trace
    let lenear_weights = Seq.map (fun (x:float*float) -> let lat,lon = x in delauney.NatNearestInterpolation(lon,lat,false,false)) grid |> Array.ofSeq
    sprintf "%A: got linear weights" sw.Elapsed |> trace
    

    let calcMean weights obs =
        Array.fold (fun state v -> let idx,w = v in state+w*obs.value.[idx] ) 0.0 weights
    
    let dist lat1 lon1 lat2 lon2= Microsoft.Research.Science.FetchClimate2.SphereMath.GetDistance(lat1,lon1,lat2,lon2)
    let sill = 227.0042+25.4346
    let cov dist = 
        let unitSph dist range =
            let tripleRange = range*range*range
            if dist>=0.0 && dist<range then
                1.5*dist/range-0.5*(dist*dist*dist)/tripleRange
            else
                1.0
        let unitNugget dist =
            if dist=0.0 then 0.0 else 1.0
        sill-(unitNugget dist)-25.4346+(unitSph dist 5263.965)*227.0042

    let calcVar targetLat targetLon weights obs = 
        let indeces = Array.map fst weights
        let w = Array.map snd weights
        GeoMath.LinearEstimator.EstimateVariance(System.Func<_,_>(cov),System.Func<_,_,_,_,_>(dist),sill,targetLat,targetLon,indeces,w,obs.lat,obs.lon)
    let means = Array.map (fun x -> calcMean x obs) lenear_weights
    sprintf "%A: mean calculated" sw.Elapsed |> trace
    let vars = Array.map2 (fun x o -> let lat,lon=o in calcVar lat lon x obs) lenear_weights grid
    sprintf "%A: variance calculated" sw.Elapsed |> trace
    sw.Stop()
    means,vars,(double(sw.ElapsedMilliseconds)*1e-03)

[<EntryPoint>]
let main argv =   
    Trace.Listeners.Add(new ConsoleTraceListener())  
    use ds = DataSet.Open "msds:csv?file=january_v3_adjusted.csv&openMode=readOnly"
    let lats  = ds.["latitude"].GetData() :?> float array
    let lons = ds.["longitude"].GetData() :?> float array
    let values = ds.["tavg_lapserate"].GetData() :?> float array
    let data = {
        lat = lats;
        lon = lons;
        value = values;
    }

    let lonsteps = [|0.25; 0.5; 0.5; 1.0; 1.0;2.0;2.0|]
    let latsteps = [| 0.25;0.25;0.5;0.5;1.0;1.0;2.0|]

    use measureDs = DataSet.Open "msds:csv?file=measurements.csv&openMode=create&appendMetadata=false"
    measureDs.IsAutocommitEnabled<-false    

    for i in 0..((Array.length lonsteps)-1) do
        let lonstep = lonsteps.[i]
        let latstep = latsteps.[i]
        let lon = [|-4.0..lonstep..36.0|]
        let lat = [|37.0..latstep..57.0|]
        let pointsCount = lon.Length*lat.Length
        let grid = [| for i in lat do yield [| for j in lon do yield (i,j)|] |] |> Seq.concat |> Seq.toArray
        let m = Array.zeroCreate subruns
        for j in 0..subruns-1 do
            let sw = Stopwatch.StartNew()
            let res = predictGrid grid data
            sw.Stop()
            let means,vars,elapsed = res
            m.[j] <- elapsed
            
            use outDs = DataSet.Open (sprintf "msds:csv?file=%d.csv&openMode=create&appendMetadata=false" pointsCount) 
            outDs.IsAutocommitEnabled <- false
            
            outDs.AddVariable<float>("mean",means,"i") |> ignore
            outDs.AddVariable<float>("var",vars,"i") |> ignore
            outDs.Commit()
        measureDs.AddVariable<float>((sprintf "%d_points" pointsCount),m,"i") |> ignore
        measureDs.Commit()
    0
    
