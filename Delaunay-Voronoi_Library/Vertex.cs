using System;
using System.Collections.Generic;

namespace Delaunay_Voronoi_Library
{
    /// <summary>
    /// The Vertex with 3D position and value
    /// </summary>
    public class Vertex
    {

        #region _fields

        private double value = double.NaN;
        private double r_max = double.NaN;
        private double r_kv = double.NaN;
        private double r_lv = double.NaN;
        private double latitude = double.NaN;
        private double longitude = double.NaN;
        private double[] Position = null;
        private List<triangle> adjacent_triangles = null;
        private VoronoiCell voronoi_cell = null;

        #endregion _fields

        #region _constructors
        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Vertex with the specified longitude, latitude and value (default equal 0).
        /// </summary>
        /// <param name="longitude">The vertex longitude.</param>
        /// <param name="latitude">The vertex latitude.</param>
        /// <param name="value">The vertex value.</param>
        public Vertex(double longitude, double latitude, double value = 0)
        {
            this.latitude = latitude;
            this.longitude = longitude;
            this.value = value;
            Position = new double[3];
            Position[0] = Math.Cos(latitude * Math.PI / 180) * Math.Cos(longitude * Math.PI / 180);
            Position[1] = Math.Cos(latitude * Math.PI / 180) * Math.Sin(longitude * Math.PI / 180);
            Position[2] = Math.Sin(latitude * Math.PI / 180);
            adjacent_triangles = new List<triangle>();
        }

        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Vertex with the specified 3D position, value (default equal 0) excluding longitude and latitude initialization.
        /// </summary>
        /// <param name="Position">The vertex position.</param>
        /// <param name="value">The vertex value .</param>
        public Vertex(double[] Position, double value = 0)
        {
            this.value = value;
            this.Position = Position;
            adjacent_triangles = new List<triangle>();
        }

        /// <summary>
        /// Initializes a new instance of the class Delaunay_Voronoi_Library.Vertex with the specified initial vertex excluding adjacent triangles.
        /// </summary>
        /// <param name="vertex">The vertex.</param>
        public Vertex(Vertex vertex)
        {
            this.latitude = vertex.latitude;
            this.longitude = vertex.longitude;
            this.value = vertex.value;
            Position = new double[3];
            Position[0] = vertex.Position[0];
            Position[1] = vertex.Position[1];
            Position[2] = vertex.Position[2];
            adjacent_triangles = new List<triangle>();
        }

        #endregion _constructors

        #region _properties
        /// <summary>
        /// Gets the X component from vertex 3D position.
        /// </summary>
        /// <value>The X component from vertex 3D position.</value>
        public double X
        {
            get
            {
                return Position[0];
            }
        }

        /// <summary>
        /// Gets the Y component from vertex 3D position.
        /// </summary>
        /// <value>The Y component from vertex 3D position.</value>
        public double Y
        {
            get
            {
                return Position[1];
            }
        }

        /// <summary>
        /// Gets the Z component from vertex 3D position.
        /// </summary>
        /// <value>The Z component from vertex 3D position.</value>
        public double Z
        {
            get
            {
                return Position[2];
            }
        }

        /// <summary>
        /// Gets the vertex latitude.
        /// </summary>
        /// <value>The vertex latitude.</value>
        public double Latitude
        {
            get
            {
                return latitude;
            }
        }

        /// <summary>
        /// Gets the vertex longitude.
        /// </summary>
        /// <value>The vertex longitude.</value>
        public double Longitude
        {
            get
            {
                return longitude;
            }
        }

        /// <summary>
        /// Gets the adjacent triangles.
        /// </summary>
        /// <value>The adjacent triangles.</value>
        public List<triangle> GetAdjacentTriangles
        {
            get
            {
                return adjacent_triangles;
            }
        }

        /// <summary>
        /// Gets or sets the vertex value.
        /// </summary>
        /// <value>The vertex value.</value>
        public double Value
        {
            get
            {
                return value;
            }
            set
            {
                this.value = value;
            }
        }

        public double GetR_MAX
        {
            get
            {
                return r_max;
            }
            set
            {
                r_max = value;
            }
        }
        public double GetR_LV
        {
            get
            {
                return r_lv;
            }
            set
            {
                r_lv = value;
            }
        }
        public double GetR_KV
        {
            get
            {
                return r_kv;
            }
            set
            {
                r_kv = value;
            }
        }
        /// <summary>
        /// Gets or sets the vertex value.
        /// </summary>
        /// <value>The vertex value.</value>
        public VoronoiCell Voronoi_Cell
        {
            get
            {
                return voronoi_cell;
            }
            set
            {
                this.voronoi_cell = value;
            }
        }
        #endregion _properties

        #region _public_methods
        /// <summary>
        /// Set new position .
        /// </summary>
        /// <param name="longitude">The vertex longitude.</param>
        /// <param name="latitude">The vertex latitude.</param>
        public void SetNewPosition(double longitude, double latitude)
        {
            this.longitude = longitude;
            this.latitude = latitude;

            Position[0] = Math.Cos(latitude * Math.PI / 180) * Math.Cos(longitude * Math.PI / 180);
            Position[1] = Math.Cos(latitude * Math.PI / 180) * Math.Sin(longitude * Math.PI / 180);
            Position[2] = Math.Sin(latitude * Math.PI / 180);

            foreach (var w in this.adjacent_triangles)
                w.Update();
        }
        #endregion _public_methods

        public Tuple<Vertex, double>[] LambdasArray = null;
    }
}
