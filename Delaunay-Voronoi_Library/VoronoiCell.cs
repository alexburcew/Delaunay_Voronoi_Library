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

        private double ca, cb, cd, cA, cB, cD, sa, sb, sd;
        private void CalcSquare()
        {
            double s = 0;
            int i = 0;
            while (i < vertices.Count - 1)
            {
                ca = (CellCentr.X * vertices[i].X + CellCentr.Y * vertices[i].Y + CellCentr.Z * vertices[i].Z);
                cb = (CellCentr.X * vertices[i + 1].X + CellCentr.Y * vertices[i + 1].Y + CellCentr.Z * vertices[i + 1].Z);
                cd = (vertices[i + 1].X * vertices[i].X + vertices[i + 1].Y * vertices[i].Y + vertices[i + 1].Z * vertices[i].Z);

                sb = Math.Sin(Math.Acos(cb));
                sa = Math.Sin(Math.Acos(ca));
                sd = Math.Sin(Math.Acos(cd));
                cA = (ca - cb * cd) / (sb * sd);
                cB = (cb - ca * cd) / (sa * sd);
                cD = (cd - cb * ca) / (sb * sa);
                s += Math.Acos(cA) + Math.Acos(cB) + Math.Acos(cD) - Math.PI;
                i++;
            }
            ca = (CellCentr.X * vertices[vertices.Count - 1].X + CellCentr.Y * vertices[vertices.Count - 1].Y + CellCentr.Z * vertices[vertices.Count - 1].Z);
            cb = (CellCentr.X * vertices[0].X + CellCentr.Y * vertices[0].Y + CellCentr.Z * vertices[0].Z);
            cd = (vertices[0].X * vertices[vertices.Count - 1].X + vertices[0].Y * vertices[vertices.Count - 1].Y + vertices[0].Z * vertices[vertices.Count - 1].Z);

            sb = Math.Sin(Math.Acos(cb));
            sa = Math.Sin(Math.Acos(ca));
            sd = Math.Sin(Math.Acos(cd));
            cA = (ca - cb * cd) / (sb * sd);
            cB = (cb - ca * cd) / (sa * sd);
            cD = (cd - cb * ca) / (sb * sa);
            s += Math.Acos(cA) + Math.Acos(cB) + Math.Acos(cD) - Math.PI;
            square = s;
        }

        #endregion _private_methods

    }
}
