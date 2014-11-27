namespace DelaunayVoronoiUnitTest

open System
open Microsoft.VisualStudio.TestTools.UnitTesting

open Delaunay_Voronoi_Library

[<TestClass>]
type UnitTest() = 
    [<TestMethod>]
    member x.TestPlanePoints () = 
        let vertices = 
            [|
                new Vertex(0.0,0.0,3.0); //testing that in first 3 points belong to one plane Delanay triangulation performed successfully
                new Vertex(90.0,0.0,6.0);
                new Vertex(180.0,0.0,9.0);
                new Vertex(270.0,0.0,12.0);
                new Vertex(0.0,90.0,15.0);
                new Vertex(0.0,-90.0,18.0)
            |]
        new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices)) |> ignore

    [<TestMethod>]
    member x.TestEqualAreasNNI () = 
        let vertices = 
            [|
                new Vertex(0.0,0.0,3.0);
                new Vertex(90.0,0.0,6.0);
                new Vertex(0.0,90.0,15.0);
                new Vertex(180.0,0.0,9.0);
                new Vertex(270.0,0.0,12.0);
                new Vertex(0.0,-90.0,18.0)
            |]
        let extractValue idx =
            vertices.[idx].Value
        let voronoi = new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices));
        let lat = Math.Asin(1.0/sqrt(3.0))*180.0/Math.PI
        let res = voronoi.NaturalNeighbourLinearCombination(45.0,lat)
        Assert.AreEqual(8.0, Array.fold (fun state t -> let idx,w = t in state+w*(extractValue idx)) 0.0 res,1e-13)

    [<TestMethod>]
    member x.TestSimpleGrid () = 
        let vertices = 
            [|
                new Vertex(5.0,40.0,1.0);
                new Vertex(30.0,40.0,2.0);
                new Vertex(5.0,50.0,3.0);
                new Vertex(30.0,50.0,4.0);                
            |]
        let extractValue idx =
            vertices.[idx].Value
        let voronoi = new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices));        
        let res = voronoi.NaturalNeighbourLinearCombination(10.0,45.0)
        let resVal = Array.fold (fun state t -> let idx,w = t in state+w*(extractValue idx)) 0.0 res        
        Assert.IsFalse (Double.IsNaN(resVal))
