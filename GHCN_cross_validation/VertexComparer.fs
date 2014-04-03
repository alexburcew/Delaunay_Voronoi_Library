module VertexComparer

open System.Collections.Generic

type Vertex = Delaunay_Voronoi_Library.Vertex

type Comparer() =
    interface IEqualityComparer<Vertex> with
        member this.Equals (v1:Vertex,v2:Vertex) =
            let treshold = 1e-5
            let res=(abs(v1.X-v2.X)<=treshold)&&(abs(v1.Y-v2.Y)<=treshold)&&(abs(v1.Z-v2.Z)<=treshold)
            res
        member this.GetHashCode (v:Vertex) =
            let hCode = int(7907.0*v.X+33550336.0*v.Y-8128.0*v.Z)
            hCode