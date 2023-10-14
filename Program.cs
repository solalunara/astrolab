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

var DKList = GetListData( "DKList" );
var DUList = GetListData( "DUList" );

for ( int i = 0; i < DUList.Length; ++i )
{
    var sr = new StreamReader( "../" + DUList[ i ] + ".txt" );
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
    ((float x, float y) LeftCut_pt, (float x, float y) RightCut_pt) = SideCreep( xydata, cutoff, interpolate: false );
    int LeftCut = xydata.ToList().IndexOf( LeftCut_pt );
    int RightCut = xydata.ToList().IndexOf( RightCut_pt );
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

    //calculate HWHM using same creep method
    float HalfMax = xydata.Max( ( (float x, float y) var ) => var.y ) / 2f;
    ((float x, float y) LeftHalf, (float x, float y) RightHalf) = SideCreep( xydata, HalfMax, interpolate: true );
    float FWHM = LeftHalf.x - RightHalf.x;
    float HWHM = FWHM / 2;

    //Shift max by background error and find new HWHM
    float HalfMaxShifted = HalfMax - SideStdDev;
    ((float x, float y) LeftHalfShift, (float x, float y) RightHalfShift) = SideCreep( xydata, HalfMaxShifted, interpolate: true );
    float FWHM_Shifted = LeftHalfShift.x - RightHalfShift.x;
    float HWHM_Shifted = FWHM_Shifted / 2;

    //Now we can use the difference as the error in the HWHM
    float HWHM_Error = MathF.Abs( HWHM_Shifted - HWHM );

    Console.WriteLine( DUList[ i ] + " - HWHM " + HWHM + " +/- " + HWHM_Error );
    sr.Close();
}

static ((float x, float y) left, (float x, float y) right) SideCreep( (float x, float y)[] xydata, float Value, bool interpolate )
{
    int left, right;
    left = right = -1;
    for ( int j = 0; j < xydata.Length / 2; ++j )
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
    (float x, float y) leftpt = (frac_l * xdiff_l + xydata[ left - 1 ].x, frac_l * ydiff_l + xydata[ left - 1 ].y);

    for ( int j = xydata.Length - 1; j >= xydata.Length / 2; --j )
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
    (float x, float y) rightpt = (frac_r * xdiff_r + xydata[ right + 1 ].x, frac_r * ydiff_r + xydata[ right + 1 ].y);

    if ( interpolate )
        return (leftpt, rightpt);
    else return (xydata[ left ], xydata[ right ]);
}