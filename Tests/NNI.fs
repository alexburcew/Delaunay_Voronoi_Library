﻿namespace DelaunayVoronoiUnitTest

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
        new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices),false) |> ignore

    [<TestMethod>]
    member x.TestEqualAreasNNI () = 
        let vertices = 
            [|
                new Vertex(0.0,0.0,3.0,0.0,0);
                new Vertex(90.0,0.0,6.0,0.0,1);
                new Vertex(0.0,90.0,15.0,0.0,2);
                new Vertex(180.0,0.0,9.0,0.0,3);
                new Vertex(270.0,0.0,12.0,0.0,4);
                new Vertex(0.0,-90.0,18.0,0.0,5)
            |]
        let extractValue idx =
            vertices.[idx].Value
        let voronoi = new Delaunay_Voronoi(new Collections.Generic.List<Vertex>(vertices),false);
        let lat = Math.Asin(1.0/sqrt(3.0))*180.0/Math.PI
        let res = voronoi.NatNearestInterpolation(45.0,lat,false,false)
        Assert.AreEqual(8.0, Array.fold (fun state t -> let idx,w = t in state+w*(extractValue idx)) 0.0 res,1e-13)
