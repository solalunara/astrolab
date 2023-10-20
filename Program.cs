#define KnownDistances
// comment/uncomment to use known vs unknown distances

static string[] GetListData( string Path )
{
    string[] List;
    using ( var List_sr = new StreamReader( Path ) )
    {
        string List_raw = List_sr.ReadToEnd();
        List_sr.BaseStream.Position = 0;
        List = new string[ List_raw.Count( ( char c ) => c == '\n' ) + 1 ];
        for ( int i = 0; i < List.Length; ++i )
        {
            List[ i ] = List_sr.ReadLine();
        }
    }
    return List;
}



#if KnownDistances
var DKList = GetListData( "DKList" );
for ( int i = 0; i < DKList.Length; ++i )
#else
var DUList = GetListData( "DUList" );
for ( int i = 0; i < DUList.Length; ++i )
#endif
{
    StreamReader sr;
    try 
    {
#if KnownDistances
        sr = new StreamReader( "../data/" + DKList[ i ] + ".txt" );
#else
        sr = new StreamReader( "../data/" + DUList[ i ] + ".txt" );
#endif
    }
    catch
    {
#if KnownDistances
        Console.WriteLine( "No data availible for " + DKList[ i ] );
#else
        Console.WriteLine( "No data availible for " + DUList[ i ] );
#endif
        continue;
    }
    string data = sr.ReadToEnd();
    string[] data_split = data.Split( new char[]{' ', '\n'}, StringSplitOptions.RemoveEmptyEntries );
    (float x,float y)[] xydata = new (float,float)[ data_split.Length / 2 ];

    for ( int j = 0; j < data_split.Length - 1; j += 2 )
        xydata[ j / 2 ] = (float.Parse( data_split[ j ] ), float.Parse( data_split[ j + 1 ] ));

    //convert frequencies to m/s
    const float c = 3e8f;
    const float l_0 = 1.42040580e9f; //rest frequency
    if ( xydata[ 0 ].x > c ) //faster than light, must be frequency
    {
        for ( int j = 0; j < xydata.Length; ++j )
            xydata[ j ].x = Math.Abs(c / l_0 * ( xydata[ j ].x - l_0 ) );
    }
    //all data now in m/s, convert all to km/s
    for ( int j = 0; j < xydata.Length; ++j )
        xydata[ j ].x = xydata[ j ].x / 1e3f;

    //separate the centre of the data from the sides
    //can't use a .where() because that might include dips
    float cutoff = xydata.Max( ( (float x, float y) var ) => var.y ) / 4f;
    ((float x, float y, int i) LeftCut_pt, (float x, float y, int i) RightCut_pt) = SideCreep( xydata, cutoff );
    int LeftCut = LeftCut_pt.i;
    int RightCut = RightCut_pt.i;
    (float, float)[] SideData_Raw = new (float, float)[ 0 ];
    if ( LeftCut != -1 )
        SideData_Raw = SideData_Raw.Union( xydata[ ..LeftCut ] ).ToArray();
    if ( RightCut != -1 )
        SideData_Raw = SideData_Raw.Union( xydata[ RightCut.. ] ).ToArray();
    float[] SideData = SideData_Raw.Select( ( (float x, float y) var ) => var.y ).ToArray();


    //use side data to calculate uncertainty in brightness
    float SideAvg = SideData.Average();
    float SideVariance = 0;
    for ( int j = 0; j < SideData.Length; ++j )
        SideVariance += MathF.Pow( SideData[ j ] - SideAvg, 2 );
    SideVariance /= SideData.Length;
    float SideStdDev = MathF.Sqrt( SideVariance );

    //seperate out all the distributions
    float HalfMax = xydata.Max( ( (float x, float y) var ) => var.y ) / 2f;
    float FWHM_Total, FWHM_Error_Total;
    FWHM_Total = FWHM_Error_Total = 0f;
    int N = 0;
    foreach ( var distribution in GetDistributions( xydata, HalfMax ) )
    {
        //calculate HWHM using same creep method
        ((float x, float y, int i) LeftHalf, (float x, float y, int i) RightHalf) = SideCreep( distribution, HalfMax );
        float FWHM = LeftHalf.x - RightHalf.x;

        //Shift max by background error and find new HWHM
        float HalfMaxShifted = HalfMax - SideStdDev;
        ((float x, float y, int i) LeftHalfShift, (float x, float y, int i) RightHalfShift) = SideCreep( distribution, HalfMaxShifted );
        float FWHM_Shifted = LeftHalfShift.x - RightHalfShift.x;

        //Now we can use the difference as the error in the FWHM
        float FWHM_Error = MathF.Abs( FWHM_Shifted - FWHM );

        FWHM_Total += FWHM;
        FWHM_Error_Total += FWHM_Error * FWHM_Error;
        ++N;
    }
    FWHM_Total /= N;
    FWHM_Error_Total /= N;
    FWHM_Error_Total = MathF.Sqrt( FWHM_Error_Total );
    if ( N > 1 )
        throw new Exception( "??? impossible state ???" );



#if KnownDistances
    Console.WriteLine( DKList[ i ] + " - FWHM " + FWHM_Total + " +/- " + FWHM_Error_Total );
#else
    Console.WriteLine( DUList[ i ] + " - FWHM " + FWHM_Total + " +/- " + FWHM_Error_Total );
#endif
    sr.Close();
}

static ((float x, float y, int i) left, (float x, float y, int i) right) SideCreep( (float x, float y)[] xydata, float Value )
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
    (float x, float y, int i) leftpt = (frac_l * xdiff_l + xydata[ left - 1 ].x, frac_l * ydiff_l + xydata[ left - 1 ].y, left);

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
    (float x, float y, int i) rightpt = (frac_r * xdiff_r + xydata[ right + 1 ].x, frac_r * ydiff_r + xydata[ right + 1 ].y, right);

    return (leftpt, rightpt);
}

static IEnumerable<(float x, float y)[]> GetDistributions( (float x, float y)[] xydata, float Value )
{
    ((float x, float y, int i) left, (float x, float y, int i) right) = SideCreep( xydata, Value );
    yield return xydata[ (left.i - 1)..(right.i + 1) ];
    /*
    for ( int i = 0; i < xydata.Length; ++i )
    {
        List<(float x, float y)> Distribution = new();
        //check for beginning of distribution
        if ( xydata[ i ].y > Value && ( ( i > 0 && xydata[ i - 1 ].y <= Value ) || i == 0 ) )
        {
            //add all of distribution above value (and ends)
            if ( i > 0 )
                Distribution.Add( xydata[ i - 1 ] );

            while ( i < xydata.Length && xydata[ i ].y > Value )
                Distribution.Add( xydata[ i++ ] );

            Distribution.Add( xydata[ i-- ] ); //set to one below so for loop inc doesn't mess us up

            //elliminate noise caused by individual points
            if ( Distribution.Count > 3 )
                yield return Distribution.ToArray();
        }
    }
    */
}