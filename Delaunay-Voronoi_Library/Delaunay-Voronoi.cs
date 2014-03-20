﻿using System;
using System.Collections.Generic;
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

        private List<Vertex> vinit = new List<Vertex>();
        public HashSet<Vertex> vertices = new HashSet<Vertex>();
        private HashSet<triangle> triangles = new HashSet<triangle>();
        private List<VoronoiCell> pol = new List<VoronoiCell>();
        private HashSet<triangle> newtriangles = new HashSet<triangle>();
        private List<triangle> pseuvotr = new List<triangle>();
        private List<Vertex> pseuvovert = new List<Vertex>();
        private HashSet<VoronoiCell> pseudopol = new HashSet<VoronoiCell>();
        private List<KeyValuePair<double, double>> apprvariogramm = null;
        #endregion _fields
        
        
        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Delaunay_Voronoi with the specified initial Delaunay_Voronoi class.
        /// </summary>
        /// <param name="delaunay_voronoi">The Delaunay_Voronoi class.</param>
        public Delaunay_Voronoi(Delaunay_Voronoi delaunay_voronoi)
        {
            Dictionary<Vertex, Vertex> dictionary_parallel_copy = new Dictionary<Vertex, Vertex>();
           
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
                
                this.triangles.Add(new triangle(vertex1,vertex2,vertex3));
            }

            foreach (var w in delaunay_voronoi.pol)
            {
                Vertex vertex;
                dictionary_parallel_copy.TryGetValue(w.GetCellCentr, out vertex);
                VoronoiCell newvc = new VoronoiCell(w,vertex);
                vertex.Voronoi_Cell = newvc;
                pol.Add(newvc);
            }
            apprvariogramm = new List<KeyValuePair<double, double>>();
            foreach(var tin in delaunay_voronoi.apprvariogramm) apprvariogramm.Add(new KeyValuePair<double,double>(tin.Key,tin.Value));
        }
        
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Delaunay_Voronoi with the specified vertices list.
        /// </summary>
        /// <param name="vertices">The vertices list.</param>
        /// <param name="Uncertainty">The uncertainty flag </param>
        public Delaunay_Voronoi(List<Vertex> vertices, bool Uncertainty = false)
        {
            int lenght = vertices.Count, i = 2;
            if (lenght < 4) { Console.WriteLine("слишком мало точек"); return; }
            double A = 0, B = 0, C = 0;
            Vertex v0 = null, v1 = null, v2 = null;
            if (lenght > 2)
            {
                v0 = vertices[0];
                v1 = vertices[1];
                do
                {
                    v2 = vertices[i++];
                    A = (v1.Y - v0.Y) * (v2.Z - v0.Z) - (v1.Z - v0.Z) * (v2.Y - v0.Y);
                    B = (v1.Z - v0.Z) * (v2.X - v0.X) - (v1.X - v0.X) * (v2.Z - v0.Z);
                    C = (v1.X - v0.X) * (v2.Y - v0.Y) - (v1.Y - v0.Y) * (v2.X - v0.X);
                } while (A == 0 && B == 0 && C == 0 && i < lenght);
            }

            if (A != 0 || B != 0 || C != 0)
            {
                var er = new double[] { v0.X + v1.X + v2.X, v0.Y + v1.Y + v2.Y, v0.Z + v1.Z + v2.Z };
                var rty = Math.Sqrt(er[0] * er[0] + er[1] * er[1] + er[2] * er[2]);
                er[0] = -er[0] / rty; er[1] =- er[1] / rty; er[2] = -er[2] / rty;
                var vini = new Vertex(er, double.NaN);

                this.vertices.Add(v0);
                this.vertices.Add(v1);
                this.vertices.Add(v2);
                this.vinit.Add(vini);
                this.vertices.Add(vini);

                triangles.Add(new triangle(v1, v2, v0));
                triangles.Add(new triangle(v0, v1, vini));
                triangles.Add(new triangle(v0, v2, vini));
                triangles.Add(new triangle(v1, v2, vini));

                foreach (var j in vertices)
                {
                    addnewpoint(j);
                }
                
                foreach (var e in vinit)
                {
                    if (TryDeletePoint(e) == false) Console.WriteLine("Невозможно удалить точку...к сожалению...");
                }

                vinit.Clear();

                Voronoi();

                var tu = DateTime.Now;
                if (Uncertainty == true) createvariogram(this.vertices.ToList(),25);
                Console.WriteLine("Variogram: {0}", DateTime.Now - tu);
            }
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
                return pol;
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
                List<Vertex> temp = new List<Vertex>();
                triangle t = q.GetAdjacentTriangles[0];
                Vertex y = t.GetVertices.First(a => a != q);
                
                for (int i = 0; i < q.GetAdjacentTriangles.Count; i++)
                {
                    temp.Add(t.GetOR);
                    y0 = y;
                    y = t.GetVertices.Single(a => (a != y) && (a != q));
                    t = q.GetAdjacentTriangles.Single(a => a.GetVertices.Contains(y) && (!a.GetVertices.Contains(y0)));
                }
                
                VoronoiCell p = new VoronoiCell(q, temp);
                pol.Add(p);
                q.Voronoi_Cell = p;
            }
        }

        void createvariogram(List<Vertex> vertices, int s)
        {
            apprvariogramm = new List<KeyValuePair<double, double>>();
            int count = vertices.Count;
            List<Stack<KeyValuePair<double, double>>> sts = new List<Stack<KeyValuePair<double, double>>>();
            List<KeyValuePair<double, double>> x = new List<KeyValuePair<double, double>>();
            
            var trt = DateTime.Now;
            Stack<KeyValuePair<double, double>> l = new Stack<KeyValuePair<double, double>>();
            
            for (int i = 0; i < count; i++)
            {
                var y = vertices[i];
                for (int j = i + 1; j < count; j++)
                {
                    var u = vertices[j];
                    double d = u.X * y.X + u.Y * y.Y + u.Z * y.Z;
                    if (d <= 0) continue;
                    l.Push(new KeyValuePair<double, double>(d, Math.Abs(u.Value - y.Value)));
                }
            }
            Console.WriteLine(DateTime.Now-trt);
            
            trt = DateTime.Now;
            x.Add(new KeyValuePair<double, double>(1, 0));
            sts.Add(l);
            
            while (sts.Count < s)
            {
                int p = sts.Max(a => a.Count);
                var q = sts.FindIndex(a => a.Count == p);
                l = sts[q];
                var xi = x[q];
                sts.Remove(l);
                x.Remove(xi);

                var xis = Math.Cos((Math.Acos(xi.Key) + Math.Acos(xi.Value)) / 2);
                x.Add(new KeyValuePair<double, double>(xi.Key, xis));
                x.Add(new KeyValuePair<double, double>(xis, xi.Value));

                var k1 = new Stack<KeyValuePair<double, double>>();
                var k2 = new Stack<KeyValuePair<double, double>>();
                sts.Add(k1); sts.Add(k2);

                for (int i = 0; i < p; i++)
                {
                    var y = l.Pop();
                    if (y.Key > xis) k1.Push(y);
                    else k2.Push(y);
                }
            }

            Console.WriteLine(DateTime.Now - trt);
            trt = DateTime.Now;
            apprvariogramm.Add(new KeyValuePair<double, double>(0, 0));
            for (int i = 0; i < s; i++)
            {
                double mo = sts[i].Sum(a => a.Value);
                mo = mo * mo / sts[i].Count;
                double d = (sts[i].Sum(v => v.Value * v.Value) - mo) / (sts[i].Count - 1);
                apprvariogramm.Add(new KeyValuePair<double, double>(Math.Acos(x[i].Value), d));
            }
            apprvariogramm.OrderByDescending(a => a.Key);
            Console.WriteLine(DateTime.Now - trt);
        }

        int addnewpoint(Vertex x, bool p = true)
        {
            foreach (var w in vinit.Where(a => Math.Abs(x.X - a.X) < 0.0000001 && Math.Abs(x.Y - a.Y) < 0.0000001 && Math.Abs(x.Z - a.Z) < 0.0000001))
            {
                w.SetNewPosition(x.Longitude, x.Latitude);
                w.Value = x.Value;
                vinit.Remove(w);
                return 0;
            }

            if (vertices.Count > 0)
            {
                var tm = vertices.ElementAt(0);
                double dis = (tm.X - x.X) * (tm.X - x.X) + (tm.Y - x.Y) * (tm.Y - x.Y) + (tm.Z - x.Z) * (tm.Z - x.Z), dis0 = 0;
                Vertex ind0 = tm;
                do
                {
                    tm = ind0;
                    foreach (var er in tm.GetAdjacentTriangles)
                    {
                        foreach (var t in er.GetVertices)
                        {
                            dis0 = (t.X - x.X) * (t.X - x.X) + (t.Y - x.Y) * (t.Y - x.Y) + (t.Z - x.Z) * (t.Z - x.Z);
                            if (dis0 < dis)
                            {
                                dis = dis0;
                                ind0 = t;
                            }
                        }
                    }
                } while (tm != ind0);

                if (dis < 0.0000000001)
                {
                    if (p == false)
                    {
                        x.Value = ind0.Value;
                    }
                    else { Console.WriteLine("Warning: SingularPoint"); }
                    return -1;
                }
            }

            triangle z = null;
            triangle l = triangles.ElementAt(0);

            if (Math.Sign(l.SignVertex(x)) == l.GetSignCenter)
            {
                l.TempIndex = x;

                double e0 = (l.GetVertices[0].X - x.X) * (l.GetVertices[0].X - x.X) + (l.GetVertices[0].Y - x.Y) * (l.GetVertices[0].Y - x.Y) + (l.GetVertices[0].Z - x.Z) * (l.GetVertices[0].Z - x.Z);
                double e1 = (l.GetVertices[1].X - x.X) * (l.GetVertices[1].X - x.X) + (l.GetVertices[1].Y - x.Y) * (l.GetVertices[1].Y - x.Y) + (l.GetVertices[1].Z - x.Z) * (l.GetVertices[1].Z - x.Z);
                double e2 = (l.GetVertices[2].X - x.X) * (l.GetVertices[2].X - x.X) + (l.GetVertices[2].Y - x.Y) * (l.GetVertices[2].Y - x.Y) + (l.GetVertices[2].Z - x.Z) * (l.GetVertices[2].Z - x.Z);

                int u = 0;
                if ((e0 >= e1) && (e2 >= e1))
                {
                    u = 1;
                }

                if ((e0 >= e2) && (e1 >= e2))
                {
                    u = 2;
                }

                while (true)
                {
                    z = null;
                    foreach (var et in l.GetVertices[u].GetAdjacentTriangles)
                    {
                        if ((et.TempIndex != x) && (Math.Sign(et.SignVertex(x)) != et.GetSignCenter))
                        {
                            z = et;
                            break;
                        }
                    }

                    if (z != null) break;

                    triangle t0 = null;
                    int u0 = 0;
                    double e = 5;
                    foreach (var et in l.GetVertices[u].GetAdjacentTriangles)
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
                    l = t0;
                    u = u0;
                }
                l = z;
            }

            HashSet<triangle> m = new HashSet<triangle>();

            anymore(m, l, x);

            if (p == false)
            {
                foreach (var d in m)
                    pseuvotr.Add(new triangle(d.GetVertices[0], d.GetVertices[1], d.GetVertices[2], false));
            }

            List<Edge> b = new List<Edge>();

            foreach (var t in m)
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
                if (p == false) newtriangles.Add(tem);
                triangles.Add(tem);
            }

            vertices.Add(x);
            return 1;
        }

        void anymore(HashSet<triangle> m, triangle l, Vertex x)
        {
            foreach (var w in l.GetVertices)
            {
                foreach (var e in w.GetAdjacentTriangles)
                {
                    if (Math.Sign(e.SignVertex(x)) != e.GetSignCenter)
                    {
                        if (!m.Contains(e))
                        {
                            m.Add(e);
                            anymore(m, e, x);
                        }
                    }
                }
            }
        }
        
        bool TryDeletePoint(Vertex v)
        {
            List<Edge> edges = new List<Edge>();
            List<Vertex> verts = new List<Vertex>();
            List<triangle> trsw = new List<triangle>();
            Dictionary<Vertex, Vertex> parallelcopy = new Dictionary<Vertex, Vertex>();
            
            foreach (var a in v.GetAdjacentTriangles)
            {
                foreach (var w in a.GetVertices)
                {
                    if (w != v && !verts.Contains(w))
                    {
                        verts.Add(w);
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
                Vertex p1=tempo[0],p2=tempo[1];
                edges.Add(new Edge(p1,p2));
            }
            vertices.Remove(v);
            Vertex v1 = verts[0],
                   v2 = edges.First(a => a.GetVertexes.Contains(v1)).GetVertexes.Single(b => b != v1),
                   v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

            int maxit = 0;
            if (verts.Count == 4) maxit = 2;
            if (verts.Count > 4) maxit = verts.Count;

            while (verts.Count > 3)
            {
                triangle t1 = new triangle(v1, v2, v3,false);
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
                        maxit--;
                        if (maxit == 0) { vertices.Add(v); return false; }
                        v1 = v2;
                        v2 = v3;
                        v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);

                    }
                    else
                    {
                        verts.Remove(v2);
                        if (verts.Count == 4) maxit = 2;
                        if (verts.Count > 4) maxit = verts.Count;
                        edges.RemoveAll(a => a.GetVertexes.Contains(v2));
                        edges.Add(new Edge(v1, v3));
                        trsw.Add(t1);
                        v2 = v3;
                        v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                    }
                }

                if (y == false)
                {
                    maxit--;
                    if (maxit == 0) { vertices.Add(v); return false; }
                    v1 = v2;
                    v2 = v3;
                    v3 = edges.Single(a => (a.GetVertexes.Contains(v2)) && (!a.GetVertexes.Contains(v1))).GetVertexes.Single(b => b != v2);
                }
            }
            var fgh = new triangle(verts[0], verts[1], verts[2],false);
            trsw.Add(fgh);

            maxit = v.GetAdjacentTriangles.Count;
            for (int k = 0; k < maxit;k++ )
            {
                var o = v.GetAdjacentTriangles[0];
                o.deltr();
                triangles.Remove(o);
            }

            foreach (var o in trsw)
            {
                triangles.Add(new triangle(o.GetVertices[0], o.GetVertices[1], o.GetVertices[2]));
            }

            return true;
        }

        void NatNearestInterpolation(Vertex vert, bool ExceptNaN, bool Uncertainty = false)
        {
            List<Vertex> lv = new List<Vertex>();
            if (apprvariogramm == null) Uncertainty = false;
            pseudopol.Clear();
            pseuvotr.Clear();
            pseuvovert.Clear();
            newtriangles.Clear();
            
            if (addnewpoint(vert, false) == -1) return;
            
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
            foreach (var q in lv)
            {
                if ((!ExceptNaN) || (ExceptNaN && !double.IsNaN(q.Value)))
                {
                    double sqss = q.Voronoi_Cell.GetSquare;
                    double sqs = pseudopol.Single(a => a.GetCellCentr == q).GetSquare;
                    double delta = sqss - sqs;
                    vert.Value += delta * q.Value;
                    sq += delta;
                    dict.Add(q, delta);
                }
            }
            vert.Value /= sq;

            if (Uncertainty == true)
            {
                double kdel = 0.0, ldel = 0.0, mdel = 0.0;
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
                    mdel = Math.Max(mdel, Math.Sqrt(((apprvariogramm[i].Value - apprvariogramm[i - 1].Value) * r + apprvariogramm[i - 1].Value * apprvariogramm[i].Key - apprvariogramm[i].Value * apprvariogramm[i - 1].Key) / (apprvariogramm[i].Key - apprvariogramm[i - 1].Key)));
                    ldel += t.Value * Math.Sqrt(((apprvariogramm[i].Value - apprvariogramm[i - 1].Value) * r + apprvariogramm[i - 1].Value * apprvariogramm[i].Key - apprvariogramm[i].Value * apprvariogramm[i - 1].Key) / (apprvariogramm[i].Key - apprvariogramm[i - 1].Key)) / sq;
                    kdel += Math.Pow(t.Value * Math.Sqrt(((apprvariogramm[i].Value - apprvariogramm[i - 1].Value) * r + apprvariogramm[i - 1].Value * apprvariogramm[i].Key - apprvariogramm[i].Value * apprvariogramm[i - 1].Key) / (apprvariogramm[i].Key - apprvariogramm[i - 1].Key)) / sq, 2);
                }
                vert.GetR_KV = Math.Sqrt(kdel);
                vert.GetR_LV = ldel;
                vert.GetR_MAX = mdel;
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
        public List<Vertex> NatNearestInterpolation(List<Vertex> f, bool parallel = false, bool ExceptNaN = true, bool Uncertainty = false)
        {
            if (parallel == false)
            {
                List<Vertex> res = new List<Vertex>();
                foreach (var w in f)
                {
                    Vertex v = new Vertex(w);
                    NatNearestInterpolation(w, ExceptNaN, Uncertainty);
                    res.Add(v);
                }
                return res;
            }
            else
            {
                int pc = Environment.ProcessorCount;

                Delaunay_Voronoi[] dvm = new Delaunay_Voronoi[pc];

                var c = DateTime.Now;
                dvm[0] = this;
                for (int i = 1; i < pc; i++)
                    dvm[i] = new Delaunay_Voronoi(this);
                Console.WriteLine("Copy {0}", DateTime.Now - c);

                int[] h = new Int32[pc + 1];
                h[0] = 0;
                h[pc] = f.Count;
                for (int i = 1; i < pc; i++)
                {
                    h[i] = (int)Math.Truncate((double)(i * f.Count / pc));
                }
                List<Vertex> res = new List<Vertex>();
                Action<int> _ForAction = (i) =>
                {
                    int _h1, _h2;
                    lock (h) { _h1 = h[i]; _h2 = h[i + 1]; }
                    List<Vertex> d = new List<Vertex>();
                    for (int j = _h1; j < _h2; j++)
                    {
                        Vertex v = new Vertex(f[j]);

                        dvm[i].NatNearestInterpolation(v, ExceptNaN, Uncertainty);

                        d.Add(v);
                    }
                    lock (res)
                    {
                        res.AddRange(d);
                    }
                };
                Parallel.For(0, pc, _ForAction);

                return res;
            }
        }

        /// <summary>
        /// Natural Nearest Interpolation.
        /// </summary>
        public List<Vertex> NatNearestInterpolation(double fromlongitude, double fromlatitude, double tolongitude, double tolatitude, int Nlongitude, int Nlatitude, bool parallel = false, bool ExceptNaN = true, bool Uncertainty = false)
        {
            List<Vertex> f = new List<Vertex>();
            for (int i = 0; i < Nlongitude; i++)
                for (int j = 0; j < Nlatitude; j++)
                    f.Add(new Vertex(fromlongitude + 1.0*(tolongitude - fromlongitude) / Nlongitude * i, fromlatitude + 1.0*(tolatitude - fromlatitude) / Nlatitude * j));
           
            return NatNearestInterpolation(f,parallel,ExceptNaN,Uncertainty);
        }
        /// <summary>
        /// Natural Nearest Interpolation.
        /// </summary>
        /// <param name="longitude">The longitude.</param>
        /// <param name="latitude">The latitude.</param>
        /// <param name="ExceptNaN">Except NaN points </param>
        /// <param name="Uncertainty">The uncertainty flag </param>
        public Vertex NatNearestInterpolation(double longitude, double latitude, bool ExceptNaN = true, bool Uncertainty = false)
        {
            Vertex w = new Vertex(longitude, latitude);
            NatNearestInterpolation(w, ExceptNaN, Uncertainty);
            return w;
        }
        #endregion _public_methods

    }
}
