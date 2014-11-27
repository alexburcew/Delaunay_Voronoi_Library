using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

namespace Delaunay_Voronoi_Library
{


    /// <summary>
    /// Delaunay_Voronoi class constructs Delaunay triangulation and Voronoi diagram.
    /// </summary>
    public class Delaunay_Voronoi
    {

        #region _fields

        public static TraceSource traceSource = new TraceSource("Delaunay_Voronoi library", SourceLevels.All);


        private List<Vertex> vinit = new List<Vertex>(); // a list of "dummy" "good positioned" NaN containing vertices to bootstrap the delaunay triangulation assembling?
        public Vertex[] Vertices;
        private HashSet<Vertex> vertices = new HashSet<Vertex>();
        private HashSet<triangle> triangles = new HashSet<triangle>();
        private List<VoronoiCell> voronoiCells = new List<VoronoiCell>();
        private HashSet<triangle> newtriangles = new HashSet<triangle>();
        private List<triangle> pseuvotr = new List<triangle>();
        private List<Vertex> pseuvovert = new List<Vertex>();
        private HashSet<VoronoiCell> pseudopol = new HashSet<VoronoiCell>();
        private List<KeyValuePair<double, double>> apprvariogramm = null;
        private List<KeyValuePair<double, double>> apprcovariogramm = null;
        #endregion _fields


        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Delaunay_Voronoi with the specified initial Delaunay_Voronoi class.
        /// </summary>
        /// <param name="delaunay_voronoi">The Delaunay_Voronoi class.</param>
        public Delaunay_Voronoi(Delaunay_Voronoi delaunay_voronoi)
        {
            Dictionary<Vertex, Vertex> dictionary_parallel_copy = new Dictionary<Vertex, Vertex>();
            int N = delaunay_voronoi.Vertices.Length;
            Vertices = new Vertex[N];
            for (int i = 0; i < N; i++)
                Vertices[i] = new Vertex(delaunay_voronoi.Vertices[i]);

            foreach (var w in delaunay_voronoi.GetVertices)
            {
                Vertex parallel_copy = new Vertex(w);
                dictionary_parallel_copy.Add(w, parallel_copy);
                this.GetVertices.Add(parallel_copy);
            }

            foreach (var w in delaunay_voronoi.triangles)
            {
                Vertex vertex1, vertex2, vertex3;

                dictionary_parallel_copy.TryGetValue(w.GetVertices[0], out vertex1);
                dictionary_parallel_copy.TryGetValue(w.GetVertices[1], out vertex2);
                dictionary_parallel_copy.TryGetValue(w.GetVertices[2], out vertex3);

                this.triangles.Add(new triangle(vertex1, vertex2, vertex3));
            }

            foreach (var w in delaunay_voronoi.voronoiCells)
            {
                Vertex vertex;
                dictionary_parallel_copy.TryGetValue(w.GetCellCentr, out vertex);
                VoronoiCell newvc = new VoronoiCell(w, vertex);
                vertex.Voronoi_Cell = newvc;
                voronoiCells.Add(newvc);
            }
            if (delaunay_voronoi.apprvariogramm != null)
            {
                apprvariogramm = new List<KeyValuePair<double, double>>();
                foreach (var tin in delaunay_voronoi.apprvariogramm) apprvariogramm.Add(new KeyValuePair<double, double>(tin.Key, tin.Value));

                apprcovariogramm = new List<KeyValuePair<double, double>>();
                foreach (var tin in delaunay_voronoi.apprcovariogramm) apprcovariogramm.Add(new KeyValuePair<double, double>(tin.Key, tin.Value));
            }

        }

        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Delaunay_Voronoi with the specified vertices list.
        /// </summary>
        /// <param name="vertices">The vertices list.</param>
        /// <param name="Uncertainty">The uncertainty flag </param>
        public Delaunay_Voronoi(List<Vertex> vertices)
        {
            int length = vertices.Count, i = 0, j = 1, k = 2;
            for (int l = 0; l < length; l++)
                vertices[l].DataIndex = l;

            if (length < 4) { traceSource.TraceEvent(TraceEventType.Critical, 1, "Too few points. Less then 4"); return; }
            Vertices = vertices.ToArray();
            double A = 0;
            Vertex v0 = null, v1 = null, v2 = null;
            while (true)
            {
                v0 = vertices[i];
                v1 = vertices[j];
                v2 = vertices[k];

                A = ((v1.Y - v0.Y) * (v2.Z - v0.Z) - (v1.Z - v0.Z) * (v2.Y - v0.Y)) * (-v0.X) +
                ((v1.X - v0.X) * (v2.Y - v0.Y) - (v1.Y - v0.Y) * (v2.X - v0.X)) * (-v0.Z) +
                 ((v1.X - v0.X) * (v2.Z - v0.Z) - (v1.Z - v0.Z) * (v2.X - v0.X)) * (-v0.Y);

                if (A != 0) break; //belong to the same great-circle??
                k++;
                if (k == length - 1) { j += 1; k = j + 1; } //this is "for loop in for loop in for loop!!!"
                if (j == length - 1) { i += 1; j = i + 1; }
                if (i == length - 1) break;

            }
            if (A != 0) //v0, v1 and v2 are not at the great circle
            {
                var er = new double[] { v0.X + v1.X + v2.X, v0.Y + v1.Y + v2.Y, v0.Z + v1.Z + v2.Z }; //vector through the triangle
                var rty = Math.Sqrt(er[0] * er[0] + er[1] * er[1] + er[2] * er[2]); //its length
                er[0] = -er[0] / rty; er[1] = -er[1] / rty; er[2] = -er[2] / rty; //normalized back vector
                var vini = new Vertex(er, double.NaN);

                this.vertices.Add(v0);
                this.vertices.Add(v1);
                this.vertices.Add(v2);
                this.vinit.Add(vini);
                this.vertices.Add(vini); //now convex hull of 4 vectors covers center of the sphere

                triangles.Add(new triangle(v1, v2, v0));
                triangles.Add(new triangle(v0, v1, vini));
                triangles.Add(new triangle(v0, v2, vini));
                triangles.Add(new triangle(v1, v2, vini));

                foreach (var jr in vertices)
                {
                    addnewpoint(jr);
                }

                foreach (var e in vinit.ToArray())
                {
                    if (TryDeletePoint(e) == false) traceSource.TraceEvent(TraceEventType.Warning, 2, "Can't remove one of the bootstrap points");
                    else
                        vinit.Remove(e);
                }

                Voronoi();

                var tu = DateTime.Now;
            }
            else
                throw new ArgumentException();
        }

        #endregion _constructors

        #region _properties
        /// <summary>
        /// Gets the voronoi cells.
        /// </summary>
        /// <value>The voronoi cells.</value>
        public List<VoronoiCell> GetVoronoiCells
        {
            get
            {
                return voronoiCells;
            }
        }
        /// <summary>
        /// Gets the triangles list.
        /// </summary>
        /// <value>The triangles.</value>
        public HashSet<triangle> GetTriangles
        {
            get
            {
                return triangles;
            }
        }
        /// <summary>
        /// Gets the vertices list.
        /// </summary>
        /// <value>The vertices list.</value>
        public HashSet<Vertex> GetVertices
        {
            get
            {
                return vertices;
            }
        }

        #endregion _properties

        #region _private_methods
        /// <summary>
        /// Constructs voronoi diagram from delaunay triangulation.
        /// </summary>
        void Voronoi()
        {

            Vertex y0;
            foreach (var q in vertices)
            {
                List<Vertex> trianglesCircumCentres = new List<Vertex>();
                triangle t = q.GetAdjacentTriangles[0]; //start from any adjecent triangle
                Vertex y = t.GetVertices.First(a => a != q); //and from any adjecent vertex

                for (int i = 0; i < q.GetAdjacentTriangles.Count; i++) //iterating over adjecent triangles, collecting their circumcentres, they form voronoi cell
                {
                    trianglesCircumCentres.Add(t.GetOR);
                    y0 = y;
                    y = t.GetVertices.Single(a => (a != y) && (a != q));
                    t = q.GetAdjacentTriangles.Single(a => a.GetVertices.Contains(y) && (!a.GetVertices.Contains(y0)));
                }

                VoronoiCell p = new VoronoiCell(q, trianglesCircumCentres);
                voronoiCells.Add(p);
                q.Voronoi_Cell = p;
            }
        }
      
        /// <summary>
        /// Adds a point to the Delaunay triangulation
        /// </summary>
        /// <param name="x"></param>
        /// <param name="isPermanentAdd">mixed responsibility here</param>
        /// <returns></returns>
        int addnewpoint(Vertex x, bool isPermanentAdd = true) 
        {
            foreach (var w in vinit.Where(a => Math.Abs(x.X - a.X) < 0.0000001 && Math.Abs(x.Y - a.Y) < 0.0000001 && Math.Abs(x.Z - a.Z) < 0.0000001)) //"dummy" bootstrap vertex hit! setting a real useful value for it. Removing from butstrap list
            {
                w.SetNewPosition(x.Longitude, x.Latitude);
                w.Value = x.Value;
                vinit.Remove(w);
                return 0;
            }

            if (vertices.Count > 0)
            {
                Vertex v0 = vertices.ElementAt(0);
                double dis = (v0.X - x.X) * (v0.X - x.X) + (v0.Y - x.Y) * (v0.Y - x.Y) + (v0.Z - x.Z) * (v0.Z - x.Z), dis0 = 0; //tracking distance to closest point
                Vertex v1 = v0;
                do
                {
                    v0 = v1;
                    foreach (var er in v0.GetAdjacentTriangles) //looking through closest point by jumping along triangles towards minimazation of distance?
                    {
                        foreach (var t in er.GetVertices)
                        {
                            dis0 = (t.X - x.X) * (t.X - x.X) + (t.Y - x.Y) * (t.Y - x.Y) + (t.Z - x.Z) * (t.Z - x.Z);
                            if (dis0 < dis)
                            {
                                dis = dis0;
                                v1 = t;
                            }
                        }
                    }
                } while (v0 != v1); //we did not move closer. The point is alread

                if (dis < 0.0000000001) //vertex hit
                {
                    if (isPermanentAdd == false) //TODO: remove mixed responsibility
                    {
                        x.Value = v1.Value;
                        x.DataIndex = v1.DataIndex;
                    }
                    else { traceSource.TraceEvent(TraceEventType.Warning, 4, "Duplicate Vertex in the initial vertex sequence. Using first one, discarding later duplicates vertex values. This can be false warning if the point is one of the initial triangle points"); }
                    return -1;
                }
            }

            triangle suspectCoveringTriangle = null;
            triangle currentIterTriangle = triangles.ElementAt(0);

            if (Math.Sign(currentIterTriangle.SignVertex(x)) == currentIterTriangle.GetSignCenter) //x is outside the triangle l, we need to find the triangle that covers vertex x
            {
                currentIterTriangle.TempIndex = x;

                //getting the direction to move in terms of l triangle vertices
                double e0 = (currentIterTriangle.GetVertices[0].X - x.X) * (currentIterTriangle.GetVertices[0].X - x.X) + (currentIterTriangle.GetVertices[0].Y - x.Y) * (currentIterTriangle.GetVertices[0].Y - x.Y) + (currentIterTriangle.GetVertices[0].Z - x.Z) * (currentIterTriangle.GetVertices[0].Z - x.Z);
                double e1 = (currentIterTriangle.GetVertices[1].X - x.X) * (currentIterTriangle.GetVertices[1].X - x.X) + (currentIterTriangle.GetVertices[1].Y - x.Y) * (currentIterTriangle.GetVertices[1].Y - x.Y) + (currentIterTriangle.GetVertices[1].Z - x.Z) * (currentIterTriangle.GetVertices[1].Z - x.Z);
                double e2 = (currentIterTriangle.GetVertices[2].X - x.X) * (currentIterTriangle.GetVertices[2].X - x.X) + (currentIterTriangle.GetVertices[2].Y - x.Y) * (currentIterTriangle.GetVertices[2].Y - x.Y) + (currentIterTriangle.GetVertices[2].Z - x.Z) * (currentIterTriangle.GetVertices[2].Z - x.Z);

                int nextVertexIndex = 0;
                if ((e0 >= e1) && (e2 >= e1))
                {
                    nextVertexIndex = 1;
                }

                if ((e0 >= e2) && (e1 >= e2))
                {
                    nextVertexIndex = 2;
                }

                while (true) //iterating through triangles, moving closer to the triangle containing x
                {
                    suspectCoveringTriangle = null;
                    foreach (var et in currentIterTriangle.GetVertices[nextVertexIndex].GetAdjacentTriangles)//checking for covering triangle among nearest neighbour trianles
                    {
                        if ((et.TempIndex != x) && (Math.Sign(et.SignVertex(x)) != et.GetSignCenter))
                        {
                            suspectCoveringTriangle = et; //triangle coveting x is found!
                            break;
                        }
                    }

                    if (suspectCoveringTriangle != null) break;

                    triangle t0 = null;
                    int u0 = 0;
                    double e = 5;
                    foreach (var et in currentIterTriangle.GetVertices[nextVertexIndex].GetAdjacentTriangles) //determining what adjecent triangle to move into next for investigation
                    {
                        if (et.TempIndex != x)
                        {
                            et.TempIndex = x;
                            e0 = (et.GetVertices[0].X - x.X) * (et.GetVertices[0].X - x.X) + (et.GetVertices[0].Y - x.Y) * (et.GetVertices[0].Y - x.Y) + (et.GetVertices[0].Z - x.Z) * (et.GetVertices[0].Z - x.Z);
                            e1 = (et.GetVertices[1].X - x.X) * (et.GetVertices[1].X - x.X) + (et.GetVertices[1].Y - x.Y) * (et.GetVertices[1].Y - x.Y) + (et.GetVertices[1].Z - x.Z) * (et.GetVertices[1].Z - x.Z);
                            e2 = (et.GetVertices[2].X - x.X) * (et.GetVertices[2].X - x.X) + (et.GetVertices[2].Y - x.Y) * (et.GetVertices[2].Y - x.Y) + (et.GetVertices[2].Z - x.Z) * (et.GetVertices[2].Z - x.Z);

                            if ((e >= e0) && (e1 >= e0) && (e2 >= e0))
                            {
                                e = e0;
                                t0 = et;
                                u0 = 0;
                            }

                            if ((e >= e1) && (e2 >= e1) && (e0 >= e1))
                            {
                                e = e1;
                                t0 = et;
                                u0 = 1;
                            }

                            if ((e >= e2) && (e0 >= e2) && (e1 >= e2))
                            {
                                e = e2;
                                t0 = et;
                                u0 = 2;
                            }
                        }
                    }
                    currentIterTriangle = t0;
                    nextVertexIndex = u0;
                }
                currentIterTriangle = suspectCoveringTriangle;
            }

            HashSet<triangle> coveringTriangles = new HashSet<triangle>(); //in case of x is at the edge of triangles???

            CheckAdjecentTrianglesForCovering(coveringTriangles, currentIterTriangle, x);

            if (isPermanentAdd == false)
            {
                foreach (var d in coveringTriangles)
                    pseuvotr.Add(new triangle(d.GetVertices[0], d.GetVertices[1], d.GetVertices[2], false));
            }

            List<Edge> b = new List<Edge>();

            foreach (var t in coveringTriangles)
            {
                if (b.RemoveAll(a => a.GetVertexes.Contains(t.GetVertices[0]) && a.GetVertexes.Contains(t.GetVertices[1])) == 0)
                {
                    b.Add(new Edge(t.GetVertices[0], t.GetVertices[1]));
                }
                if (b.RemoveAll(a => a.GetVertexes.Contains(t.GetVertices[2]) && a.GetVertexes.Contains(t.GetVertices[1])) == 0)
                {
                    b.Add(new Edge(t.GetVertices[2], t.GetVertices[1]));
                }
                if (b.RemoveAll(a => a.GetVertexes.Contains(t.GetVertices[2]) && a.GetVertexes.Contains(t.GetVertices[0])) == 0)
                {
                    b.Add(new Edge(t.GetVertices[2], t.GetVertices[0]));
                }
                t.deltr();
                triangles.Remove(t);
            }
            foreach (var q in b)
            {
                var tem = new triangle(x, q.GetVertexes[0], q.GetVertexes[1]);
                if (isPermanentAdd == false) newtriangles.Add(tem);
                triangles.Add(tem);
            }

            vertices.Add(x);
            return 1;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="coveringTriangles"></param>
        /// <param name="startTriangle"></param>
        /// <param name="x"></param>
        void CheckAdjecentTrianglesForCovering(HashSet<triangle> coveringTriangles, triangle startTriangle, Vertex x)
        {
            foreach (var w in startTriangle.GetVertices)
            {
                foreach (var e in w.GetAdjacentTriangles)
                {
                    if (Math.Sign(e.SignVertex(x)) != e.GetSignCenter)
                    {
                        if (!coveringTriangles.Contains(e))
                        {
                            coveringTriangles.Add(e);
                            CheckAdjecentTrianglesForCovering(coveringTriangles, e, x);
                        }
                    }
                }
            }
        }

        static double InnerProduct(Vertex v1,Vertex v2)
        {
            return v1.X * v2.X + v1.Y * v2.Y + v1.Z * v2.Z;
        }

        bool TryDeletePoint(Vertex v)
        {
            List<Edge> unenclosedEdges = new List<Edge>();
            List<Vertex> unenclosedVerts = new List<Vertex>();
            List<triangle> newlyConstructedTriangles = new List<triangle>();
            Dictionary<Vertex, Vertex> parallelcopy = new Dictionary<Vertex, Vertex>();

            bool isConvexHullVertex = vertices.Any(v5 => v5 != v && InnerProduct(v,v5)<=0.0);
            if (isConvexHullVertex) //removingf the vertex will cause the convex hull not to cover 0,0,0 point
            {
                traceSource.TraceEvent(TraceEventType.Error,5,"Unable to remove the vertex as this causes triangulation not to cover 0,0,0 point");
                return false;
            }


            foreach (var a in v.GetAdjacentTriangles)
            {
                foreach (var w in a.GetVertices)
                {
                    if (w != v && !unenclosedVerts.Contains(w))
                    {
                        unenclosedVerts.Add(w);
                    }
                }
                List<Vertex> tempo = new List<Vertex>();
                foreach (var w in a.GetVertices)
                {
                    if (w != v)
                    {
                        tempo.Add(w);
                    }
                }
                Vertex p1 = tempo[0], p2 = tempo[1];
                unenclosedEdges.Add(new Edge(p1, p2));
            }
            vertices.Remove(v);
            Vertex v1 = unenclosedVerts[0],
                   v2 = unenclosedEdges.First(a => a.GetVertexes.Contains(v1)).GetVertexes.Single(b => b != v1),
                   v3 = unenclosedEdges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

            int maxit = 0;
            if (unenclosedVerts.Count == 4) maxit = 2;
            if (unenclosedVerts.Count > 4) maxit = unenclosedVerts.Count;

            while (unenclosedVerts.Count > 3)
            {
                triangle t1 = new triangle(v1, v2, v3, false); //triangle under investigation
                bool y = true;
                foreach (var v4 in vertices)
                {
                    int c = Math.Sign(t1.SignVertex(v4));
                    if ((c != 0) && (c != t1.GetSignCenter))
                    {
                        y = false;
                        break;
                    }
                }

                if (y == true) //good. all vertices are at the same side from the t1 plane. it is probably dalaynay triangle
                {
                    var g = new triangle(v1, v3, new Vertex(new double[] { 0, 0, 0 }), false);
                    if (Math.Sign(g.SignVertex(v2)) == Math.Sign(g.SignVertex(v))) //check for delaunay properies?
                    {
                        maxit--;
                        if (maxit == 0) { vertices.Add(v); return false; }
                        v1 = v2;
                        v2 = v3;
                        v3 = unenclosedEdges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

                    }
                    else
                    {
                        unenclosedVerts.Remove(v2);
                        if (unenclosedVerts.Count == 4) maxit = 2;
                        if (unenclosedVerts.Count > 4) maxit = unenclosedVerts.Count;
                        unenclosedEdges.RemoveAll(a => a.GetVertexes.Contains(v2));
                        unenclosedEdges.Add(new Edge(v1, v3));
                        newlyConstructedTriangles.Add(t1);
                        v2 = v3;
                        v3 = unenclosedEdges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                    }
                }
                else //poor t1 combination it s covered by convex hull
                {
                    maxit--;
                    if (maxit == 0) { vertices.Add(v); return false; }
                    v1 = v2;
                    v2 = v3;
                    v3 = unenclosedEdges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                }
            }
            var fgh = new triangle(unenclosedVerts[0], unenclosedVerts[1], unenclosedVerts[2], false);
            newlyConstructedTriangles.Add(fgh);

            maxit = v.GetAdjacentTriangles.Count;
            for (int k = 0; k < maxit; k++)
            {
                var o = v.GetAdjacentTriangles[0];
                o.deltr();
                triangles.Remove(o);
            }

            foreach (var o in newlyConstructedTriangles)
            {
                triangles.Add(new triangle(o.GetVertices[0], o.GetVertices[1], o.GetVertices[2]));
            }

            return true;
        }

        void DeletePoint3(Vertex v, bool p)
        {
            List<Edge> edges = new List<Edge>();
            List<Vertex> verts = new List<Vertex>();
            List<triangle> triangls = new List<triangle>();

            vertices.Remove(v);
            foreach (var a in v.GetAdjacentTriangles)
            {
                triangls.Add(a);
                a.GetVertices.Remove(v);
                foreach (var q0 in a.GetVertices) q0.GetAdjacentTriangles.Remove(a);
                edges.Add(new Edge(a.GetVertices[0], a.GetVertices[1]));
                foreach (var l in a.GetVertices) if (!verts.Contains(l)) verts.Add(l);
            }

            v.GetAdjacentTriangles.Clear();

            pseuvovert = new List<Vertex>(verts);
            foreach (var r in triangls)
                triangles.Remove(r);

            triangls.Clear();

            Vertex v1 = verts[0],
                   v2 = edges.First(a => a.GetVertexes.Contains(v1)).GetVertexes.Single(b => b != v1),
                   v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

            while (verts.Count > 3)
            {
                triangle t1 = new triangle(v1, v2, v3, true);
                bool y = true;
                foreach (var v4 in vertices)
                {
                    int c = Math.Sign(t1.SignVertex(v4));
                    if ((c != 0) && (c != t1.GetSignCenter))
                    {
                        y = false;
                        break;
                    }
                }

                if (y == true)
                {
                    var g = new triangle(v1, v3, new Vertex(new double[] { 0, 0, 0 }), false);
                    if (Math.Sign(g.SignVertex(v2)) == Math.Sign(g.SignVertex(v)))
                    {
                        t1.deltr();
                        v1 = v2;
                        v2 = v3;
                        v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

                    }
                    else
                    {
                        verts.Remove(v2);
                        edges.RemoveAll(a => a.GetVertexes.Contains(v2));
                        edges.Add(new Edge(v1, v3));
                        triangles.Add(t1);
                        v2 = v3;
                        v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                    }
                }

                if (y == false)
                {
                    t1.deltr();
                    v1 = v2;
                    v2 = v3;
                    v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                }
            }
            var fgh = new triangle(verts[0], verts[1], verts[2], true);
            triangles.Add(fgh);
        }

        void crossdel(List<Vertex> vx)
        {
            pseudopol.Clear();
            foreach (var q in vx)
            {

                triangle t = q.GetAdjacentTriangles[0];
                Vertex y = t.GetVertices.First(a => a != q);
                Vertex y0;
                List<Vertex> j = new List<Vertex>();
                for (int i = 0; i < q.GetAdjacentTriangles.Count; i++)
                {
                    j.Add(t.GetOR);
                    y0 = y;
                    y = t.GetVertices.Single(a => (a != y) && (a != q));
                    t = q.GetAdjacentTriangles.Single(a => a.GetVertices.Contains(y) && (!a.GetVertices.Contains(y0)));
                }
                VoronoiCell p = new VoronoiCell(q, j);
                pseudopol.Add(p);
            }
        }
        Tuple<int, double>[] NatNearestInterpolation(Vertex vert, bool CheckForNaN)
        {
            bool Uncertainty = false;
            List<Vertex> lv = new List<Vertex>();            
            pseudopol.Clear();
            pseuvotr.Clear();
            pseuvovert.Clear();
            newtriangles.Clear();

            if (addnewpoint(vert, false) == -1) //exact match to one of the initial nodes
            {
                return new Tuple<int, double>[] { Tuple.Create(vert.DataIndex, 1.0) };
            }

            foreach (var e1 in pseuvotr)
                foreach (var e2 in e1.GetVertices)
                {
                    if (!lv.Contains(e2)) lv.Add(e2);
                }
            lv.Add(vert);

            foreach (var q in lv)
            {
                List<Vertex> temp = new List<Vertex>();
                triangle t = q.GetAdjacentTriangles[0];
                Vertex y = t.GetVertices.First(a => a != q);
                Vertex y0;
                for (int i = 0; i < q.GetAdjacentTriangles.Count; i++)
                {
                    temp.Add(t.GetOR);
                    y0 = y;
                    y = t.GetVertices.Single(a => (a != y) && (a != q));
                    t = q.GetAdjacentTriangles.Single(a => a.GetVertices.Contains(y) && (!a.GetVertices.Contains(y0)));
                }
                VoronoiCell p = new VoronoiCell(q, temp);
                pseudopol.Add(p);
            }
            lv.Remove(vert);

            double sq = 0;
            Dictionary<Vertex, double> dict = new Dictionary<Vertex, double>();
            Dictionary<Vertex, double> deltas = new Dictionary<Vertex, double>();
            foreach (var q in lv)
            {
                if ((!CheckForNaN) || (CheckForNaN && !double.IsNaN(q.Value)))
                {
                    double sqss = q.Voronoi_Cell.GetSquare;
                    double sqs = pseudopol.Single(a => a.GetCellCentr == q).GetSquare;
                    double delta = sqss - sqs;
                    deltas.Add(q, delta);
                    vert.Value += delta * q.Value;
                    sq += delta;
                    dict.Add(q, delta);
                }
            }
            vert.Value /= sq;

            if (Uncertainty == true)
            {
                double kt = 0.0, lt = 0.0, mt = 0.0;
                foreach (var t in dict)
                {
                    var r = Math.Acos(t.Key.X * vert.X + t.Key.Y * vert.Y + t.Key.Z * vert.Z);
                    var x = 0.0;
                    int i = -1;
                    do
                    {
                        i++;
                        x = apprvariogramm[i].Key;
                    }
                    while (x < r);
                    kt = ((apprcovariogramm[i].Value - apprcovariogramm[i - 1].Value) * r + apprcovariogramm[i - 1].Value * apprcovariogramm[i].Key - apprcovariogramm[i].Value * apprcovariogramm[i - 1].Key) / (apprcovariogramm[i].Key - apprcovariogramm[i - 1].Key);
                    lt += t.Value / sq * t.Key.Sigma * (kt + Math.Sqrt(kt * kt + Math.Sqrt(((apprvariogramm[i].Value - apprvariogramm[i - 1].Value) * r + apprvariogramm[i - 1].Value * apprvariogramm[i].Key - apprvariogramm[i].Value * apprvariogramm[i - 1].Key) / (apprvariogramm[i].Key - apprvariogramm[i - 1].Key)) / t.Key.Sigma / t.Key.Sigma - 1));
                }

            }

            foreach (var n in newtriangles)
            {
                n.deltr();
                triangles.Remove(n);
            }
            foreach (var n in pseuvotr)
            {
                var tem = new triangle(n.GetVertices[0], n.GetVertices[1], n.GetVertices[2]);
                triangles.Add(tem);
            }
            vertices.Remove(vert);

            return deltas.Select(entry => Tuple.Create(entry.Key.DataIndex, entry.Value / sq)).ToArray();
        }

        #endregion _private_methods

        #region _public_methods
        /// <summary>
        /// Natural Nearest Interpolation.
        /// </summary>
        /// <param name="verts">The vertices.</param>
        /// <param name="parallel">parallel flag.</param>
        /// <param name="ExceptNaN">Except NaN points </param>
        /// <param name="Uncertainty">The uncertainty flag </param>
        public List<Vertex> NatNearestInterpolation(List<Vertex> f,bool checkForNaN)
        {            

                List<Vertex> res = new List<Vertex>();
                foreach (var w in f)
                {
                    Vertex v = new Vertex(w);
                    NatNearestInterpolation(w, checkForNaN);
                    res.Add(v);
                }
                return res;
            
        }

        /// <summary>
        /// Natural Nearest Interpolation.
        /// </summary>
        public List<Vertex> NaturalNeighbourInterpolation(double fromlongitude, double fromlatitude, double tolongitude, double tolatitude, int Nlongitude, int Nlatitude, bool checkforNaN = true)
        {
            bool parallel = false,  Uncertainty = false;
            List<Vertex> f = new List<Vertex>();
            for (int i = 0; i < Nlongitude; i++)
                for (int j = 0; j < Nlatitude; j++)
                    f.Add(new Vertex(fromlongitude + 1.0 * (tolongitude - fromlongitude) / Nlongitude * i, fromlatitude + 1.0 * (tolatitude - fromlatitude) / Nlatitude * j));

            return NatNearestInterpolation(f, checkforNaN);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="longitude">The longitude.</param>
        /// <param name="latitude">The latitude.</param>
        /// <param name="ExceptNaN">Except NaN points </param>
        /// <param name="Uncertainty">The uncertainty flag </param>
        public Tuple<int, double>[] NaturalNeighbourLinearCombination(double longitude, double latitude, bool checkforNaN = true)
        {
            bool ExceptNaN = true, Uncertainty = false;
            Vertex w = new Vertex(longitude, latitude);
            var res = NatNearestInterpolation(w, checkforNaN);
            return res;
        }

        #endregion _public_methods

    }
}
