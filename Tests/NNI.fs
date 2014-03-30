namespace UnitTestProject1

open System
open Microsoft.VisualStudio.TestTools.UnitTesting

open Delaunay_Voronoi_Library

[<TestClass>]
type UnitTest() = 
    [<TestMethod>]
    member x.TestEqualAreasNNI () = 
        let vertices = 
            [|
                new Vertex(0.0,0.0,3.0);
                new Vertex(90.0,0.0,6.0);
                new Vertex(180.0,0.0,9.0);
                new Vertex(270.0,0.0,12.0);
                new Vertex(0.0,90.0,15.0);
                new Vertex(0.0,-90.0,18.0)
            |]
        let voronoi = new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices),false);
        let res = voronoi.NatNearestInterpolation(45.0,45.0,false,false)
        Assert.AreEqual(8.0, res.Value,1e-12)
