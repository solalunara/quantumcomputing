using System.Reflection.Metadata.Ecma335;

namespace quantumbits;

using System.Globalization;
using System.Numerics;
using System.Runtime.Intrinsics.X86;
using static complex_utils;

public struct qbit
{
    internal qbit( complex alpha, complex beta )
    {
        this.alpha = alpha;
        this.beta = beta;
    }

    internal complex alpha;
    internal complex beta;

    public readonly double Probability_0() => Sqr( alpha );
    public readonly double Probability_1() => Sqr( beta );

    public bool Collapse()
    {
        Random rnd = new();
        if ( rnd.NextDouble() < Probability_1() )
        {
            beta = 1;
            alpha = 0;
            return true;
        }
        else
        {
            beta = 0;
            alpha = 1;
            return false;
        }
    }

    public override string? ToString()
    {
        return alpha + " |0> + " + beta + " |1>";
    }

    public void LeftRotation()
    {
        complex alpha = ( this.alpha - this.beta ) / MathF.Sqrt( 2 );
        complex beta = ( this.alpha + this.beta ) / MathF.Sqrt( 2 );
        this.alpha = alpha;
        this.beta = beta;
    }
    public void RightRotation()
    {
        complex alpha = ( this.alpha + this.beta ) / MathF.Sqrt( 2 );
        complex beta = (-this.alpha + this.beta ) / MathF.Sqrt( 2 );
        this.alpha = alpha;
        this.beta = beta;
    }
    public void S_Shift()
    {
        complex alpha = i * this.alpha;
        this.alpha = alpha;
    }
    public void T_Shift()
    {
        complex alpha = -this.alpha;
        complex beta = -i * this.beta;
        this.alpha = alpha;
        this.beta = beta;
    }

    public static implicit operator qbit( bool i )
    {
        return i switch
        {
            false => new(1, 0),
            true => new(0, 1),
        };
    }
    public static qbit_group operator ^( qbit a, qbit b )
    {
        complex alpha = a.alpha * b.alpha;
        complex beta = a.alpha * b.beta;
        complex gamma = a.beta * b.alpha;
        complex delta = a.beta * b.beta;

        (gamma, delta) = (delta, gamma);

        return new qbit_group( alpha, beta, gamma, delta );
    }
}

public struct qbit_group
{
    internal qbit_group( params complex[] a  )
    {
        int n = a.Length;
        while ( ( n >>= 1 ) >= 1 ) N++;
        this.a = a;
    }
    internal qbit_group( params qbit[] q )
    {
        N = q.Length;
        a = NestedArrayMultip( q ).ToArray();
    }
    readonly int N = 0;
    internal complex[] a;

    private IEnumerable<complex> NestedArrayMultip( qbit[] q )
    {
        if ( q.Length == 1 )
        {
            yield return q[ 0 ].alpha;
            yield return q[ 0 ].beta;
        }
        else 
        {
            foreach ( var z in NestedArrayMultip( q[ 1.. ] ) )
                yield return q[ 0 ].alpha * z;
            foreach ( var z in NestedArrayMultip( q[ 1.. ] ) )
                yield return q[ 0 ].beta * z;
        }
    }
    private IEnumerable<bool> IntToBools( int i )
    {
            for ( int j = 0; j < N; ++j )
                yield return ( i & ( 1 << ( j ) ) ) != 0;
    }

    //converts symbols (e.g. |101>) to indexes (e.g. 5 for |101>)
    public complex GetVal( int i ) => GetVal( IntToBools( i ).ToArray() );
    public readonly complex GetVal( params bool[] i ) => a[ GetValIndex( i ) ];
    public int GetValIndex( int i ) => GetValIndex( IntToBools( i ).ToArray() );
    public readonly int GetValIndex( params bool[] i )
    {
        if ( i.Length != N )
            throw new InvalidOperationException( "wrong number of symbols (" + i.Length + ") for group size (" + N + ")" );

        int n = 0;
        for ( int j = 0; j < N; ++j )
            if( i[ j ] ) 
                n |= 1 << j;
        return n;
    }

    public override string? ToString()
    {
        string? s = "";
        uint n = 0;
        for ( int i = 0; i < a.Length; ++i )
        {
            if ( i != 0 )
                s += " + ";

            char[] symbols = new char[ N ];
            for ( int j = 0; j < N; ++j )
                symbols[ N - 1 - j ] = ( n & ( 1 << j ) ) == 0 ? '0' : '1';
            ++n;
            
            s += a[ i ] + " |" + new string( symbols ) + ">";
        }
        return s;
    }

    public qbit_group Xor( int ControlBit, int TargetBit )
    {
        if ( TargetBit == ControlBit )
            throw new InvalidOperationException( "Target bit and control bit cannot be the same" );

        matrix xor_1 = matrix.Identity( 2 );
        matrix xor_2 = TargetBit == 0 ? matrix.X() : matrix.Identity( 2 );
        if ( ControlBit == 0 )
        {
            xor_1 = matrix.ZeroZero();
            xor_2 = matrix.OneOne();
        }
        for ( int i = 1; i < N; ++i )
        {
            if ( i == ControlBit )
            {
                xor_1 %= matrix.ZeroZero();
                xor_2 %= matrix.OneOne();
            }
            else if ( i == TargetBit )
            {
                xor_1 %= matrix.Identity( 2 );
                xor_2 %= matrix.X();
            }
            else
            {
                xor_1 %= matrix.Identity( 2 );
                xor_2 %= matrix.Identity( 2 );
            }
        }
        matrix xor = xor_1 + xor_2;

        return xor * this;
    }
    public qbit_group RightShift( int n )
    {
        matrix m;
        if ( n == 0 )
            m = matrix.RightRotation2x2();
        else
            m = matrix.Identity( 2 );
        for ( int i = 1; i < N; ++i )
            if ( i == n ) m %= matrix.RightRotation2x2();
            else m %= matrix.Identity( 2 );
        return m * this;
    }
    public qbit_group LeftShift( int n )
    {
        matrix m;
        if ( n == 0 )
            m = matrix.LeftRotation2x2();
        else
            m = matrix.Identity( 2 );
        for ( int i = 1; i < N; ++i )
            if ( i == n ) m %= matrix.LeftRotation2x2();
            else m %= matrix.Identity( 2 );
        return m * this;
    }
    public qbit_group SShift( int n )
    {
        matrix m;
        if ( n == 0 )
            m = matrix.SShift2x2();
        else
            m = matrix.Identity( 2 );
        for ( int i = 1; i < N; ++i )
            if ( i == n ) m %= matrix.SShift2x2();
            else m %= matrix.Identity( 2 );
        return m * this;
    }
    public qbit_group TShift( int n )
    {
        matrix m;
        if ( n == 0 )
            m = matrix.TShift2x2();
        else
            m = matrix.Identity( 2 );
        for ( int i = 1; i < N; ++i )
            if ( i == n ) m %= matrix.TShift2x2();
            else m %= matrix.Identity( 2 );
        return m * this;
    }

    public void TryDisentangle( out qbit?[] qbits )
    {
        qbits = new qbit?[ N ];
        for ( int i = N; --i >= 0; )
        {
            bool separable = true;
            int n_0 = 0;
            int n_1 = 1 << i;
            complex frac = -double.MaxValue;
            complex a = 0;
            complex b = 0;
            while ( n_0 < 1 << i && separable )
            {
                complex v_0 = GetVal( n_0 );
                complex v_1 = GetVal( n_1 );

                a += Sqr( v_0 );
                b += Sqr( v_1 );

                if ( v_0 != 0 )
                {
                    complex new_frac = v_0 / v_1;
                    if ( frac != -double.MaxValue )
                        if ( new_frac != frac ) separable = false;
                    frac = v_0 / v_1;
                }

                n_1++;
                n_0++;
            }
            double arg = Math.Atan2( frac.b, frac.a );
            b = new( b.a * Math.Cos( arg ), b.a * Math.Sin( arg ) );
            if ( separable )
            {
                //(complex a, complex b) = EstimateAB( n.ToArray() );
                qbits[ N - 1 - i ] = new( a, b );
            }
            else
                qbits[ N - 1 - i ] = null;
        }
    }

}

public struct matrix
{
    internal matrix( complex[,] data )
    {
        this.data = data;
    }
    private complex[,] data;

    public complex this[ int i, int j ]
    {
        get { return data[ i, j ]; }
        set { data[ i, j ] = value; }
    }
    public static implicit operator matrix( qbit_group q )
    {
        complex[,] data = new complex[ q.a.Length, 1 ];
        for ( int i = 0; i < q.a.Length; ++i )
            data[ i, 0 ] = q.a[ i ];
        return new matrix( data ).Round();
    }
    public static implicit operator qbit_group( matrix m )
    {
        complex[] data = new complex[ m.data.GetLength( 0 ) ];
        for ( int i = 0; i < data.Length; ++i )
            data[ i ] = m[ i, 0 ];
        return new qbit_group( data );
    }
    public static matrix operator *( matrix m1, matrix m2 )
    {
        if ( m1.data.GetLength( 1 ) != m2.data.GetLength( 0 ) )
            throw new InvalidOperationException( "Matrix dimensions are not compatible for multiplication." );

        int rowsA = m1.data.GetLength( 0 );
        int colsA = m1.data.GetLength( 1 );
        int colsB = m2.data.GetLength( 1 );

        complex[,] result = new complex[ rowsA, colsB ];

        for ( int i = 0; i < rowsA; i++ )
        {
            for ( int j = 0; j < colsB; j++ )
            {
                complex sum = 0;
                for ( int k = 0; k < colsA; k++ )
                    sum += m1[ i, k ] * m2[ k, j ];
                result[ i, j ] = sum;
            }
        }
        return new matrix( result ).Round();
    }
    public static matrix operator +( matrix m1, matrix m2 )
    {
        if ( m1.data.GetLength( 0 ) != m2.data.GetLength( 0 ) || 
            m1.data.GetLength( 1 ) != m2.data.GetLength( 1 ) )
            throw new InvalidOperationException( "Matrix dimensions are not compatible for multiplication." );

        complex[,] result = new complex[ m1.data.GetLength( 0 ), m1.data.GetLength( 1 ) ];
        for ( int i = 0; i < m1.data.GetLength( 0 ); ++i )
            for ( int j = 0; j < m1.data.GetLength( 1 ); ++j )
                result[ i, j ] = m1[ i, j ] + m2[ i, j ];
        return new matrix( result ).Round();
    }

    //tensor product
    public static matrix operator %( matrix m1, matrix m2 )
    {
        int rowsA = m1.data.GetLength( 0 );
        int colsA = m1.data.GetLength( 1 );
        int rowsB = m2.data.GetLength( 0 );
        int colsB = m2.data.GetLength( 1 );

        complex[,] result = new complex[ rowsA * rowsB, colsA * colsB ];

        for ( int i = 0; i < rowsA; i++ )
            for ( int j = 0; j < colsA; j++ )
                for ( int k = 0; k < rowsB; k++ )
                    for ( int l = 0; l < colsB; l++ )
                        result[ i * rowsB + k, j * colsB + l ] = m1[ i, j ] * m2[ k, l ];

        return new matrix( result ).Round();
    }

    public matrix Round()
    {
        for ( int i = 0; i < data.GetLength( 0 ); ++i )
            for ( int j = 0; j < data.GetLength( 1 ); ++j )
                data[ i, j ] = new complex( Math.Round( data[ i, j ].a, 5 ), Math.Round( data[ i, j ].b, 5 ) );
        return this;
    }

    public static matrix Identity( int n )
    {
        complex[,] data = new complex[ n, n ];
        for ( int i = 0; i < n; ++i )
            data[ i, i ] = 1;
        return new( data );
    }
    public static matrix LeftRotation2x2()
    {
        double s = 1 / Math.Sqrt( 2 );
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = s;
        data[ 1, 0 ] = s;
        data[ 0, 1 ] = -s;
        data[ 1, 1 ] = s;
        return new( data );
    }
    public static matrix RightRotation2x2()
    {
        double s = 1 / Math.Sqrt( 2 );
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = s;
        data[ 1, 0 ] = -s;
        data[ 0, 1 ] = s;
        data[ 1, 1 ] = s;
        return new( data );
    }
    public static matrix SShift2x2()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = i;
        data[ 1, 0 ] = 0;
        data[ 0, 1 ] = 0;
        data[ 1, 1 ] = 1;
        return new( data );
    }
    public static matrix TShift2x2()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = -1;
        data[ 1, 0 ] = 0;
        data[ 0, 1 ] = 0;
        data[ 1, 1 ] = -i;
        return new( data );
    }
    public static matrix X()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = 0;
        data[ 1, 0 ] = 1;
        data[ 0, 1 ] = 1;
        data[ 1, 1 ] = 0;
        return new( data );
    }
    public static matrix Y()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = 0;
        data[ 1, 0 ] = i;
        data[ 0, 1 ] = -i;
        data[ 1, 1 ] = 0;
        return new( data );
    }
    public static matrix Z()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = 1;
        data[ 1, 0 ] = 0;
        data[ 0, 1 ] = 0;
        data[ 1, 1 ] = -1;
        return new( data );
    }
    public static matrix ZeroZero()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = 1;
        data[ 1, 0 ] = 0;
        data[ 0, 1 ] = 0;
        data[ 1, 1 ] = 0;
        return new( data );
    }
    public static matrix OneOne()
    {
        complex[,] data = new complex[ 2, 2 ];
        data[ 0, 0 ] = 0;
        data[ 1, 0 ] = 0;
        data[ 0, 1 ] = 0;
        data[ 1, 1 ] = 1;
        return new( data );
    }
}

public struct complex : IEquatable<complex>, IFormattable
{
    internal complex( double a, double b )
    {
        this.a = a;
        this.b = b;
    }

    internal double a = 0;
    internal double b = 0;

    public override string? ToString()
    {
        if ( b == 0 )
            return a.ToString();
        if ( a == 0 )
            return b + "i";
        return a.ToString() + " + " + b + "i";
    }

    public bool Equals(complex other) => this == other;

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        return this.ToString() ?? "";
    }

    public static implicit operator complex( double a ) => new( a, 0 );

    public static complex operator +( complex z ) => z;
    public static complex operator -( complex z ) => new( -z.a, -z.b );
    public static complex operator ~( complex z ) => new( z.a, -z.b );
    public static complex operator +( complex z1, complex z2 ) => new( z1.a + z2.a, z1.b + z2.b );
    public static complex operator -( complex z1, complex z2 ) => new( z1.a - z2.a, z1.b - z2.b );
    public static complex operator *( complex z1, complex z2 ) => new( z1.a * z2.a - z1.b * z2.b, z1.a * z2.b + z1.b * z2.a );
    public static complex operator /( complex z, double f ) => new( z.a / f, z.b / f );
    public static complex operator /( complex z1, complex z2 ) 
    {
        if ( z2.b == 0 && z1.b == 0 )
            return new( z1.a / z2.a, 0 );
        else if ( z2.a == 0 && z1.a == 0 )
            return new( 0, z1.b / z2.b );
        else
            return new( ( z1 * ~z2 ).a / ( z2 * ~z2 ).a, ( z1 * ~z2 ).b / ( z2 * ~z2 ).b );
    }

    public static bool operator ==( complex z1, complex z2 ) => z1.a == z2.a && z1.b == z2.b;
    public static bool operator !=( complex z1, complex z2 ) => z1.a != z2.a || z1.b != z2.b;

    public static implicit operator bool( complex z ) => z != 0;
}

public static class complex_utils
{
    public static double Re( complex z ) => z.a;
    public static double Im( complex z ) => z.b;
    public static double Sqr( complex z ) => ( z * ~z ).a;
    public static complex i = new( 0, 1 );
}