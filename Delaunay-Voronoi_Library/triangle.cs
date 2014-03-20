using System;
using System.Collections.Generic;

namespace Delaunay_Voronoi_Library
{
    /// <summary>
    /// The triangle 
    /// </summary>
    public class triangle
    {

        #region _fields

        private List<Vertex> vertices = null;
        private int signcenter;
        private Vertex OR = null;
        private double A, B, C;
        private Vertex tempindex = null;

        #endregion _fields

        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.triangle with the specified tops and connection option (default equal true).
        /// </summary>
        /// <param name="vertex1">The first vertex.</param>
        /// <param name="vertex2">The second vertex.</param>
        /// <param name="vertex3">The third vertex.</param>
        /// <param name="connection">The connection option.</param>
        public triangle(Vertex vertex1, Vertex vertex2, Vertex vertex3, bool connection = true)
        {
            vertices = new List<Vertex>();

            vertices.Add(vertex1);
            vertices.Add(vertex2);
            vertices.Add(vertex3);

            if (connection == true)
            {
                vertex1.GetAdjacentTriangles.Add(this);
                vertex2.GetAdjacentTriangles.Add(this);
                vertex3.GetAdjacentTriangles.Add(this);
            }

            Calculate();

        }
        #endregion _constructors

        #region _properties
        /// <summary>
        /// Gets the triangle vertices list.
        /// </summary>
        /// <value>The triangle vertices list.</value>
        public List<Vertex> GetVertices
        {
            get
            {
                return vertices;
            }
        }

        /// <summary>
        /// Gets the sign of sphere center relative to the plane of the triangle.
        /// </summary>
        /// <value>The sign of sphere center relative to the plane of the triangle.</value>
        public int GetSignCenter
        {
            get
            {
                return signcenter;
            }
        }

        /// <summary>
        /// Gets the circumcenter of the triangle.
        /// </summary>
        /// <value>The circumcenter of the triangle.</value>
        public Vertex GetOR
        {
            get
            {
                if (OR == null) CalcOR();
                return OR;
            }
        }

        /// <summary>
        /// Gets and sets the temporary index of this triangle.
        /// </summary>
        /// <value>The temporary index of this triangle.</value>
        public Vertex TempIndex
        {
            get
            {
                return tempindex;
            }
            set
            {
                tempindex = value;
            }
        }
        #endregion _properties

        #region _private_methods

        private void Calculate()
        {
            A = (vertices[1].Y - vertices[0].Y) * (vertices[2].Z - vertices[0].Z) - (vertices[1].Z - vertices[0].Z) * (vertices[2].Y - vertices[0].Y);
            B = (vertices[1].Z - vertices[0].Z) * (vertices[2].X - vertices[0].X) - (vertices[1].X - vertices[0].X) * (vertices[2].Z - vertices[0].Z);
            C = (vertices[1].X - vertices[0].X) * (vertices[2].Y - vertices[0].Y) - (vertices[1].Y - vertices[0].Y) * (vertices[2].X - vertices[0].X);

            signcenter = Math.Sign(SignVertex(new Vertex(new double[] { 0, 0, 0 })));
        }

        private void CalcOR()
        {

            double lc01 = Math.Sqrt((vertices[0].X - vertices[1].X) * (vertices[0].X - vertices[1].X) +
                                    (vertices[0].Y - vertices[1].Y) * (vertices[0].Y - vertices[1].Y) +
                                    (vertices[0].Z - vertices[1].Z) * (vertices[0].Z - vertices[1].Z));

            double la21 = Math.Sqrt((vertices[2].X - vertices[1].X) * (vertices[2].X - vertices[1].X) +
                                    (vertices[2].Y - vertices[1].Y) * (vertices[2].Y - vertices[1].Y) +
                                    (vertices[2].Z - vertices[1].Z) * (vertices[2].Z - vertices[1].Z));

            double lb02 = Math.Sqrt((vertices[0].X - vertices[2].X) * (vertices[0].X - vertices[2].X) +
                                    (vertices[0].Y - vertices[2].Y) * (vertices[0].Y - vertices[2].Y) +
                                    (vertices[0].Z - vertices[2].Z) * (vertices[0].Z - vertices[2].Z));

            double polup = (lc01 + lb02 + la21) / 2;
            double sq2 = polup * (polup - lc01) * (polup - lb02) * (polup - la21);

            double alfaA = la21 * la21 / 8 / sq2
                          * ((vertices[0].X - vertices[1].X) * (vertices[0].X - vertices[2].X) +
                             (vertices[0].Y - vertices[1].Y) * (vertices[0].Y - vertices[2].Y) +
                             (vertices[0].Z - vertices[1].Z) * (vertices[0].Z - vertices[2].Z));

            double alfaB = lb02 * lb02 / 8 / sq2
                          * ((vertices[1].X - vertices[0].X) * (vertices[1].X - vertices[2].X) +
                             (vertices[1].Y - vertices[0].Y) * (vertices[1].Y - vertices[2].Y) +
                             (vertices[1].Z - vertices[0].Z) * (vertices[1].Z - vertices[2].Z));

            double alfaC = lc01 * lc01 / 8 / sq2
                          * ((vertices[2].X - vertices[0].X) * (vertices[2].X - vertices[1].X) +
                             (vertices[2].Y - vertices[0].Y) * (vertices[2].Y - vertices[1].Y) +
                             (vertices[2].Z - vertices[0].Z) * (vertices[2].Z - vertices[1].Z));

            double Ox = alfaA * vertices[0].X + alfaB * vertices[1].X + alfaC * vertices[2].X;
            double Oy = alfaA * vertices[0].Y + alfaB * vertices[1].Y + alfaC * vertices[2].Y;
            double Oz = alfaA * vertices[0].Z + alfaB * vertices[1].Z + alfaC * vertices[2].Z;

            double k = Math.Sqrt(Ox * Ox + Oy * Oy + Oz * Oz);
            OR = new Vertex(new double[] { Ox / k, Oy / k, Oz / k });

        }
        
        #endregion _private_methods

        #region _public_methods
        /// <summary>
        /// Gets the sign of vertex relative to the plane of the triangle.
        /// </summary>
        /// <param name="vertex">The vertex.</param>
        public int SignVertex(Vertex vertex)
        {
            if ((vertices[0] == vertex) || (vertices[1] == vertex) || (vertices[2] == vertex))
                return 0;

            return Math.Sign((vertex.X - vertices[0].X) * A +
                   (vertex.Y - vertices[0].Y) * B +
                   (vertex.Z - vertices[0].Z) * C);
        }

        /// <summary>
        /// Update the triangle parametrs.
        /// </summary>
        public void Update()
        {
            Calculate();

            if (OR != null)
                CalcOR();
        }

        /// <summary>
        /// Disconnect this triangle.
        /// </summary>
        public void deltr()
        {
            foreach (var d in vertices)
                d.GetAdjacentTriangles.Remove(this);
        }

        #endregion _public_methods
    }
}
