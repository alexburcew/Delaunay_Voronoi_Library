using System;
using System.Collections.Generic;

namespace Delaunay_Voronoi_Library
{
    /// <summary>
    /// The Voronoi Cell with the centr and vertexes list 
    /// </summary>
    public class VoronoiCell
    {
     
        #region _fields

        private double square = -1;
        private Vertex CellCentr = null;
        private List<Vertex> vertices = null;
        
        #endregion _fields

        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.VoronoiCell with the specified initial VoronoiCell and new central vertex.
        /// </summary>
        /// <param name="Voronoi_Cell">The Voronoi cell.</param>
        /// <param name="centr">The Voronoi cell central vertex.</param>
        public VoronoiCell(VoronoiCell Voronoi_Cell, Vertex centr)
        {
            this.CellCentr = centr;
            this.square = Voronoi_Cell.square;

            this.vertices = new List<Vertex>();
            foreach (var w in Voronoi_Cell.vertices)
                this.vertices.Add(new Vertex(w));
        }

        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.VoronoiCell with the specified cell center and the cell vertices list.
        /// </summary>
        /// <param name="CellCentr">The central vertex.</param>
        /// <param name="vertices">The cell vertices list.</param>
        public VoronoiCell(Vertex CellCentr, List<Vertex> vertices)
        {
            this.CellCentr = CellCentr;
            this.vertices = vertices;
        }

        #endregion _constructors

        #region _properties
        /// <summary>
        /// Gets the cell center.
        /// </summary>
        /// <value>The cell center.</value>
        public Vertex GetCellCentr
        {
            get
            {
                return CellCentr;
            }
        }
        /// <summary>
        /// Gets vertices.
        /// </summary>
        /// <value>The vertices.</value>
        public List<Vertex> GetVertices
        {
            get
            {
                return vertices;
            }
        }
        /// <summary>
        /// Gets the cell square.
        /// </summary>
        /// <value>The cell square.</value>
        public double GetSquare
        {
            get
            {
                if (square == -1) CalcSquare();
                return square;
            }
        }
        #endregion _properties

        #region _private_methods

        private void CalcSquare()
        {
            double s = 0.0;
            for (int i = 0; i < vertices.Count - 1; i++) // composing total area of triangles sum. Fan of triangles with the common vertex - the centre of the cell
            {
                Vertex v1 = vertices[i];
                Vertex v2 = vertices[i + 1];
                Vertex v3 = CellCentr;

                //defining 2 vectors - the 2 sides of triangles
                double bx = v1.X - v3.X, by = v1.Y - v3.Y, bz = v1.Z - v3.Z;
                double cx = v2.X - v3.X, cy = v2.Y - v3.Y, cz = v2.Z - v3.Z;

                //cross product vector
                // A = B x C
                //see http://en.wikipedia.org/wiki/Cross_product
                double ax = by * cz - bz * cy, ay = bz * cx - bx * cz, az = bx * cy - by * cx;

                double parallelogramArea = Math.Sqrt(ax * ax + ay * ay + az * az);
                s += parallelogramArea;
            }
            square = s * 0.5;
        }
        #endregion _private_methods

    }
}
