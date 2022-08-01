const RIVERS = Dict(
    :Skellefte => ([
        # name, approximate capacity, reservoir size, reported head, reservoir level limits (high & low)
        # Due to confidential data, the numbers on delay been replaced with 2h for all plants. 
        #                    MW    HE     m      m      m
        Plant(:Rebnis,       64, 205560,  82.0, 513.0, 499.5), # Upstream plant:
        Plant(:Sädva,        31, 168050,  45.0, 477.0, 460.7), # name, dischargedelay (h), spilldelay (h), flowshare (default=1.0)
        Plant(:Hornavan,      0, 209170,     0, 426.1, 423.5,  Upstream(:Rebnis,       2, 2),
                                                               Upstream(:Sädva,        2, 2)),
        Plant(:Bergnäs,       8, 216110,   4.7, 420.1, 418.0,  Upstream(:Hornavan,     2, 2)),
        Plant(:Slagnäs,       7,   0.03,   5.0, 414.4, 413.8,  Upstream(:Bergnäs,      2, 2)), 
        Plant(:Bastusel,    100,   8200,  71.5, 408.5, 407.5,  Upstream(:Slagnäs,      2, 2)),
        Plant(:Grytfors,     31,   1250,  21.8, 332.0, 331.0,  Upstream(:Bastusel,     2, 2)),
        Plant(:Gallejaur,   214,   3500,  78.8, 310.0, 307.5,  Upstream(:Grytfors,     2, 2)),
        Plant(:Vargfors,    131,   4000,  49.2, 230.5, 228.5,  Upstream(:Gallejaur,    2, 2)),
        Plant(:Rengård,      36,   1300,  19.3, 181.0, 180.0,  Upstream(:Vargfors,     2, 2)),
        Plant(:Båtfors,      44,   1330,  16.9, 161.5, 160.5,  Upstream(:Rengård,      2, 2)),
        Plant(:Finnfors,     44,    300,  20.7, 144.2, 143.2,  Upstream(:Båtfors,      2, 2)),
        Plant(:Granfors,     40,    280,  18.6, 123.5, 122.5,  Upstream(:Finnfors,     2, 2)), 
        Plant(:Krångfors,    62,    330,  29.9, 104.7, 103.8,  Upstream(:Granfors,     2, 2)), 
        Plant(:Selsfors,     61,    500,  22.2,  74.6,  72.6,  Upstream(:Krångfors,    2, 2)),
        Plant(:Kvistforsen, 130,   1120,  50.6,  52.0,  50.5,  Upstream(:Selsfors,     2, 2)), 
        Plant(:Bergsby,       0,    180,     0,     1,   0.1,  Upstream(:Kvistforsen,  2, 2)) 
    ],


    Dict( #Number of turbines per plant
        :Rebnis => [1],
        :Sädva => [1],
        :Hornavan => [],
        :Bergnäs => [1],
        :Slagnäs => [1],
        :Bastusel => [1],
        :Grytfors => [1],
        :Gallejaur => [1,2],
        :Vargfors => [1,2],
        :Rengård => [1],
        :Båtfors => [1,2],
        :Finnfors => [1,2],
        :Granfors => [1,2],
        :Krångfors => [1,2,3],
        :Selsfors => [1,2],
        :Kvistforsen => [1,2],
        :Bergsby => []
    ),

    Dict( #Data points on etafunction for each turbine. d=discharge, e=eta. 
    #Due to confidential data, the numbers on the efficiency points have been replaced with the same efficiency points for all plants. 
        (:Rebnis,1)     =>  [(d=58.0, e=0.85),(d=67.5, e=0.92),(d=85.0, e=0.84)],
        (:Sädva,1)      =>  [(d=55.0, e=0.85),(d=63.0, e=0.92),(d=81.5, e=0.84)],
        (:Bergnäs,1)    =>  [(d=73.0, e=0.85),(d=100.0, e=0.92),(d=125.0, e=0.84)],
        (:Slagnäs,1)    =>  [(d=70.0, e=0.85),(d=100.0, e=0.92),(d=120.0, e=0.84)],
        (:Bastusel,1)   =>  [(d=95.0, e=0.85),(d=120.0, e=0.92),(d=140.0, e=0.84)],
        (:Grytfors,1)   =>  [(d=56.0, e=0.85),(d=120.0, e=0.92),(d=175.0, e=0.84)],
        (:Gallejaur,1)  =>  [(d=90.0, e=0.85),(d=114.0, e=0.92),(d=140.0, e=0.84)],
        (:Gallejaur,2)  =>  [(d=105.0, e=0.85),(d=155.0, e=0.92),(d=165.0, e=0.84)],
        (:Vargfors,1)   =>  [(d=68.0, e=0.85),(d=110.0, e=0.92),(d=132.0, e=0.84)],
        (:Vargfors,2)   =>  [(d=69.0, e=0.85),(d=117.0, e=0.92),(d=164.0, e=0.84)],
        (:Rengård,1)    =>  [(d=55.0, e=0.85),(d=100.0, e=0.92),(d=224.0, e=0.84)],
        (:Båtfors,1)    =>  [(d=35.0, e=0.85),(d=98.0, e=0.92),(d=165.0, e=0.84)],
        (:Båtfors,2)    =>  [(d=35.0, e=0.85),(d=97.0, e=0.92),(d=165.0, e=0.84)],
        (:Finnfors,1)   =>  [(d=32.0, e=0.85),(d=70.0, e=0.92),(d=133.0, e=0.84)],
        (:Finnfors,2)   =>  [(d=32.0, e=0.85),(d=60.0, e=0.92),(d=90.0, e=0.84)],
        (:Granfors,1)   =>  [(d=35.0, e=0.85),(d=68.0, e=0.92),(d=124.0, e=0.84)],
        (:Granfors,2)   =>  [(d=37.0, e=0.85),(d=60.0, e=0.92),(d=85.0, e=0.84)],
        (:Krångfors,1)  =>  [(d=22.5, e=0.85),(d=40.0, e=0.92),(d=62.0, e=0.84)],
        (:Krångfors,2)  =>  [(d=40.0, e=0.85),(d=76.0, e=0.92),(d=113.5, e=0.84)],
        (:Krångfors,3)  =>  [(d=50.0, e=0.85),(d=56.0, e=0.92),(d=62.0, e=0.84)],
        (:Selsfors,1)   =>  [(d=31.0, e=0.85),(d=75.0, e=0.92),(d=125.0, e=0.84)],
        (:Selsfors,2)   =>  [(d=69.0, e=0.85),(d=124.0, e=0.92),(d=192.0, e=0.84)],
        (:Kvistforsen,1) => [(d=37.0, e=0.85),(d=120.0, e=0.92),(d=162.0, e=0.84)],
        (:Kvistforsen,2) => [(d=37.0, e=0.85),(d=120.0, e=0.92),(d=162.0, e=0.84)],
    ),

    function (year, minspill, minlevel, minflow)
        # Due to confidential data, the numbers on minspill and reservoir levels have been replaced with arbitrary numbers, and only a few plants are left as examples while the rest are removed.
        minspill[:, :Sädva] .= 1

        minspill[:, :Hornavan] .= 10
        minspill[timerange(year, "05-01", "05-31"), :Hornavan] .= 10
        minspill[timerange(year, "06-01", "09-15"), :Hornavan] .= 10
        minspill[timerange(year, "09-16", "09-30"), :Hornavan] .= 10
        minlevel[timerange(year, "06-23", "08-15"), :Hornavan, :forebay] .= 426
    
        minspill[:, :Bergnäs] .= 10    
        minspill[timerange(year, "05-15", "05-31"), :Bergnäs] .= 10
        minspill[timerange(year, "06-01", "06-14"), :Bergnäs] .= 10
        minspill[timerange(year, "06-15", "07-31"), :Bergnäs] .= 10
        minspill[timerange(year, "08-01", "08-14"), :Bergnäs] .= 10
        minspill[timerange(year, "08-15", "08-31"), :Bergnäs] .= 10
    
        minflow[:, :Kvistforsen] .= 20
        minspill[:, :Bergsby] .= 10
    end)
)

function riverdata(rivername)
    if rivername in keys(RIVERS)
        return RIVERS[rivername] # Tuple: (plants vector, set minlimits function)
    else
        error("No data for $rivername river network.")
    end
end
