﻿//#define KnownDistances
// comment to read python output, uncomment to generate python input

using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;

static IEnumerable<string> GetListData( string Path )
{
    using var List_sr = new StreamReader(Path);
    while ( List_sr.Peek() != -1 )
        yield return List_sr.ReadLine() ?? throw new ArgumentNullException(nameof(Path), "File has an empty line");
}

static Dictionary<string, statval> GetStatVarData( string Path )
{
    StreamReader sr = new( Path );
    string data = sr.ReadToEnd();
    string[] data_split = data.Split( new char[]{' ', '\n', ',', '\t'}, StringSplitOptions.RemoveEmptyEntries );
    (string name, statval val)[] incl_data = new (string, statval)[ data_split.Length / 3 ];
    for ( int i = 0; i < data_split.Length; i += 3 )
    {
        incl_data[ i / 3 ].name = data_split[ i ];
        incl_data[ i / 3 ].val = new statval( float.Parse( data_split[ i + 1 ] ), float.Parse( data_split[ i + 2 ] ) );
    }

    return incl_data.ToDictionary( ((string name, statval incl) a) => a.name, ((string name, statval incl) a) => new statval(a.incl.val, a.incl.err) );
}

static IEnumerable<(pt[] x, string name)> GetXYDataForEachFile( IEnumerable<string> list )
{
    foreach ( string file in list )
    {
        StreamReader sr;
        try 
        {
            sr = new StreamReader( "data/" + file + ".txt" );
        }
        catch
        {
            Console.WriteLine( "No data availible for " + file );
            continue;
        }
        string data = sr.ReadToEnd();
        string[] data_split = data.Split( new char[]{' ', '\n', ',', '\t'}, StringSplitOptions.RemoveEmptyEntries );
        (pt[] x, string name) xydata;
        xydata.x = new pt[ data_split.Length / 2 ];
        xydata.name = file;

        for ( int j = 0; j < data_split.Length - 1; j += 2 )
            xydata.x[ j / 2 ] = (float.Parse( data_split[ j ] ), float.Parse( data_split[ j + 1 ] ));

        //convert frequencies to m/s
        const float c = 3e8f;
        const float l_0 = 1.42040580e9f; //rest frequency
        if ( xydata.x[ 0 ].x > c ) //faster than light, must be frequency
        {
            for ( int j = 0; j < xydata.x.Length; ++j )
                xydata.x[ j ].x = Math.Abs(c / l_0 * ( xydata.x[ j ].x - l_0 ) );
        }
        //all data now in m/s, convert all to km/s
        for ( int j = 0; j < xydata.x.Length; ++j )
            xydata.x[ j ].x = xydata.x[ j ].x / 1e3f;

        yield return xydata;
    }
}

static float GetStdDev( pt[] xydata )
{
    //separate the centre of the data from the sides
    //can't use a .where() because that might include dips
    float cutoff = xydata.Max( ( pt var ) => var.y ) / 4f;
    ((pt x, int i) LeftCut_pt, (pt x, int i) RightCut_pt) = SideCreep( xydata, cutoff );
    int LeftCut = LeftCut_pt.i;
    int RightCut = RightCut_pt.i;
    pt[] SideData_Raw = new pt[ 0 ];
    if ( LeftCut != -1 )
        SideData_Raw = SideData_Raw.Union( xydata[ ..LeftCut ] ).ToArray();
    if ( RightCut != -1 )
        SideData_Raw = SideData_Raw.Union( xydata[ RightCut.. ] ).ToArray();
    float[] SideData = SideData_Raw.Select( ( pt var ) => var.y ).ToArray();


    //use side data to calculate uncertainty in brightness
    float SideAvg = SideData.Average();
    float SideVariance = 0;
    for ( int j = 0; j < SideData.Length; ++j )
        SideVariance += MathF.Pow( SideData[ j ] - SideAvg, 2 );
    SideVariance /= SideData.Length;
    return MathF.Sqrt( SideVariance );
}

static statval GetHWHM( pt[] xydata, float SideStdDev )
{
    //seperate out all the distributions
    float HalfMax = xydata.Max( ( pt var ) => var.y ) / 2f;
    var distribution = GetDistribution( xydata, HalfMax );

    //calculate HWHM using same creep method
    ((pt x, int i) LeftHalf, (pt x, int i) RightHalf) = SideCreep( distribution, HalfMax );
    float HWHM = ( LeftHalf.x.x - RightHalf.x.x ) / 2;

    //Shift max by background error and find new HWHM
    float HalfMaxShifted = HalfMax - SideStdDev;
    ((pt x, int i) LeftHalfShift, (pt x, int i) RightHalfShift) = SideCreep( distribution, HalfMaxShifted );
    float HWHM_Shifted = ( LeftHalfShift.x.x - RightHalfShift.x.x ) / 2;

    //Now we can use the difference as the error in the HWHM
    float HWHM_Error = MathF.Abs( HWHM_Shifted - HWHM );

    return (HWHM, HWHM_Error);
}

static statval GetWeightedMean( pt[] xydata )
{
    float num = 0;
    float dnm = 0;
    for ( int i = 0; i < xydata.Length; ++i )
    {
        num += xydata[ i ].x * xydata[ i ].y;
        dnm += xydata[ i ].y;
    }
    return (num/dnm, 1/dnm);
}

static ((pt x, int i) left, (pt x, int i) right) SideCreep( pt[] xydata, float Value )
{
    int left, right;
    left = right = -1;
    for ( int j = 0; j < xydata.Length; ++j )
    {
        if ( xydata[ j ].y > Value )
        {
            left = j;
            break;
        }
    }
    if ( left == 0 ) left = 1;
    float ydiff_l = xydata[ left ].y - xydata[ left - 1 ].y;
    float xdiff_l = xydata[ left - 1 ].x - xydata[ left ].x;
    float part_ydiff_l = Value - xydata[ left - 1 ].y;
    float frac_l = part_ydiff_l / ydiff_l;
    frac_l = (frac_l < -1) ? -1 : ((frac_l > 1) ? 1 : frac_l); //clamp
    (pt x, int i) leftpt = ((frac_l * xdiff_l + xydata[ left - 1 ].x, frac_l * ydiff_l + xydata[ left - 1 ].y), left);

    for ( int j = xydata.Length - 1; j >= 0; --j )
    {
        if ( xydata[ j ].y > Value )
        {
            right = j;
            break;
        }
    }
    if ( right == xydata.Length - 1 ) right = xydata.Length - 2;
    float ydiff_r = xydata[ right ].y - xydata[ right + 1 ].y;
    float xdiff_r = xydata[ right ].x - xydata[ right + 1 ].x;
    float part_ydiff_r = Value - xydata[ right + 1 ].y;
    float frac_r = part_ydiff_r / ydiff_r;
    frac_r = (frac_r < -1) ? -1 : ((frac_r > 1) ? 1 : frac_r); //clamp
    (pt x, int i) rightpt = ((frac_r * xdiff_r + xydata[ right + 1 ].x, frac_r * ydiff_r + xydata[ right + 1 ].y), right);

    return (leftpt, rightpt);
}

static pt[] GetDistribution( pt[] xydata, float Value )
{
    ((pt x, int i) left, (pt x, int i) right) = SideCreep( xydata, Value );
    return xydata[ (left.i - 1)..(right.i + 1) ];
}

static statval GetRotSpeed( statval HWHM, statval incl )
{
    //convert to radians
    incl.val = incl.val * MathF.PI / 180f;
    incl.err = incl.err * MathF.PI / 180f;

    float sin = MathF.Sin( incl.val );    
    float sin_uncert = MathF.Cos( incl.val ) * incl.err;

    float RotSpeed = HWHM.val / sin;
    float RotSpeed_Error = RotSpeed * MathF.Sqrt( MathF.Pow( HWHM.err / HWHM.val, 2 ) + MathF.Pow( sin_uncert / sin, 2 ) );
    return (RotSpeed, RotSpeed_Error);
}

#if KnownDistances
static statval GetDistance( statval distance_modulus )
{
    double dist = Math.Pow( 10, 1 + distance_modulus.val / 5d );
    double dist_err = distance_modulus.err * Math.Log( 10d ) / 5d * dist;

    //pc to m
    dist *= 3.0857e16d;
    dist_err *= 3.0857e16d;

    return ((float)dist, (float)dist_err);
}
#endif

/// MAIN CODE ///
var inclinations = GetStatVarData( "inclinations" );
var fluxdensities = GetStatVarData( "fluxdensities" );
var distance_moduli = GetStatVarData( "distance-moduli" );

//data manip
#if KnownDistances
StreamWriter Outfile = new( "data.csv" );
Outfile.WriteLine( "Data - RotVelocity - RotVelocity Error - Luminosity - Luminosity Error - Distance - Distance Error - Velocity - Velocity Error" );
foreach ( (pt[] x, string name) in GetXYDataForEachFile( GetListData( "DKList" ) ) )
#else
StreamWriter Outfile = new( "data_predicted.csv" );
Outfile.WriteLine( "Data - RotVelocity - RotVelocity Error - Luminosity - Luminosity Error - Distance - Distance Error - Velocity - Velocity Error" );
//statval A = (4.217957e+01f, 7.783958e+00f);
//statval B = (4.050670926540902f, 1.4944752485181898f);
statval B = (7.735729998874423f, 0.4130204204558944f);
statval A = (6.0506603477263695e-18f, 1.3208230713361719e-17f);
foreach ( (pt[] x, string name) in GetXYDataForEachFile( GetListData( "DUList" ) ) )
#endif
{
    statval HWHM = GetHWHM( x, GetStdDev( x ) );
    statval WeightedMean = GetWeightedMean( x );
    statval incl = inclinations[ name ];
    statval RotSpeed = GetRotSpeed( HWHM, incl );
    statval FluxDensity = fluxdensities[ name ];

#if KnownDistances
    statval dist = GetDistance( distance_moduli[ name ] );
    double distf = dist.val * 3.2408E-23;
    double dist_errf = dist.err * 3.2408E-23;
    double lum = FluxDensity.val * 4d * Math.PI * dist.val * dist.val;
    double lum_err = lum * Math.Sqrt( Math.Pow( FluxDensity.err / FluxDensity.val, 2 ) + Math.Pow( dist.err / dist.val, 2 ) );

    lum /= 1e27;
    lum_err /= 1e27;

    Outfile.WriteLine( name + "," + RotSpeed.val + "," + RotSpeed.err + "," + lum + "," + lum_err  + "," + distf + "," + dist_errf + "," + WeightedMean.val + "," + WeightedMean.err );

#else
    double lum = A.val * Math.Pow( RotSpeed.val, B.val );
    double lum_err = Math.Pow( A.err / A.val, 2 ) + Math.Pow( RotSpeed.err * B.val / RotSpeed.val, 2 ) + Math.Pow( B.err * Math.Log( RotSpeed.val ), 2 );
    lum_err = lum * Math.Sqrt( lum_err );

    lum *= 1e27;
    lum_err *= 1e27;

    //double lum_log = A.val + B.val * Math.Log( RotSpeed.val );
    //double lum_log_low = ( A.val - A.err ) + ( B.val + B.err ) * Math.Log( RotSpeed.val - RotSpeed.err );
    //double lum_log_high = ( A.val + A.err ) + ( B.val - B.err ) * Math.Log( RotSpeed.val + RotSpeed.err );
    //double lum_log_err = 0.5 * ( lum_log_high - lum_log_low );
    
    //double lum = Math.Exp( lum_log );
    //double lum_err = lum * lum_log_err;

    double dist = Math.Sqrt( lum / ( 4 * Math.PI * FluxDensity.val ) );
    double dist_err = Math.Pow( .5f * dist * lum_err / lum, 2 ) + Math.Pow( .5f * dist * FluxDensity.err / FluxDensity.val, 2 );
    dist_err = Math.Sqrt( dist_err );

    dist *= 3.2408E-23;
    dist_err *= 3.2408E-23;
    
    Outfile.WriteLine( name + "," + RotSpeed.val + "," + RotSpeed.err + "," + lum * 1e-27 + "," + lum_err * 1e-27 + "," + dist + "," + dist_err + "," + WeightedMean.val + "," + WeightedMean.err );

    dist_err /= 2;

    Console.WriteLine( name + ": " + dist + " +/- " + dist_err );
#endif
}
Outfile.Close();


// structures
struct pt
{
    public pt( float x, float y )
    {
        this.x = x;
        this.y = y;
    }
    public float x;
    public float y;

    public static implicit operator (float x, float y)( pt a ) => ( a.x, a.y );
    public static implicit operator pt( (float x, float y) p ) => new( p.x, p.y );
}
struct statval
{
    public statval( float val, float err )
    {
        this.val = val;
        this.err = err;
    }
    public float val;
    public float err;

    public static implicit operator (float val, float err)( statval a ) => ( a.val, a.err );
    public static implicit operator statval( (float val, float err) p ) => new( p.val, p.err );
}