using XLSX, JLD2, Dates

function readmatrix(file, sheet, range; swap=nothing)
    mat = XLSX.readdata(file, sheet, range)
    swap !== nothing && replace!(mat, swap)
    return replace(float.(mat), missing => -Inf)
end
readmatrix(file, sheet, lastcol, year) = readmatrix(file, sheet, "B2:$lastcol$(daysinyear(year)*24 + 1)")

function read_raw_data()
    println("Reading elspot prices...")
    pricefolder = "$DATAFOLDER/Elspot_prices"
    elspot_price = Array{Float64}(undef, 0, 1)
    for year = 2004:2010        #OBS fixa daylightsaving hours. Double the price for one hour, and one hour has 0 in price
        println(year)
        dr = Date(year,1,1):Day(1):Date(year,12,31)
        file = "stosek"*last("$year", 2)*".xlsx"
        price_data = XLSX.readdata("$pricefolder/1996-2010/$file", "stosek"*last("$year", 2), "A1:Y400")
        price_data[ismissing.(price_data)].=0
        check_date=[price_data[i,1] in dr for i in 1:400]
        price_data = price_data[check_date,:]
        price_data=reshape(price_data[:,2:end],length(price_data[:,2:end]),1)
        elspot_price=vcat(elspot_price, price_data)
    end
    for year = 2011:2012        #OBS fixa daylightsaving hours. Double the price for one hour, and one hour has 0 in price
        println(year)
        dr = Date(year,1,1):Day(1):Date(year,12,31)
        file = "lulsek"*last("$year", 2)*".xlsx"
        price_data = XLSX.readdata("$pricefolder/2011-2012/$file", "lulsek"*last("$year", 2), "A1:Y400")
        price_data[ismissing.(price_data)].=0
        check_date=[price_data[i,1] in dr for i in 1:400]
        price_data = price_data[check_date,:]
        price_data=reshape(price_data[:,2:end],length(price_data[:,2:end]),1)
        elspot_price=vcat(elspot_price, price_data)
    end
    for year = 2013:2020
        println(year)
        nr = daysinyear(year)*24+4
        file = "elspot-prices_$(year)_hourly_sek.xlsx"
        price_data = XLSX.readdata("$pricefolder/2013-2020/$file", "elspot-prices_$(year)_hourly_sek", "D4:D$nr")
        price_data = vec(price_data)
        deleteat!(price_data, findall(ismissing, price_data))
        elspot_price=vcat(elspot_price, price_data)
    end

    folder = "$DATAFOLDER/Rådata"
    println("Reading Skelleftekraft production data...")
    prod = Array{Float64}(undef, 0, 11)
    for year = 2004:2020
        println(year)
        mat = readmatrix("$folder/Prod_skeälven_$year.xlsx", "Blad1", "M", year)
        prod = vcat(prod, mat[:, [1; 2; 4:12]])
    end

    println("Reading Skelleftekraft water level data...")
    hourly = Array{Float64}(undef, 0, 38)
    for year = 2004:2020
        println(year)
        mat = readmatrix("$folder/Timdata_skeälven_$year.xlsx", "Blad1", "AM", year)
        hourly = vcat(hourly, mat)
    end
    level_upper, level_lower = hourly[:,[15;17;21:2:end]], hourly[:,[16;19;22:2:end]]

    println("Reading Skelleftekraft discharge and spill data...")
    skk_discharge = Array{Float64}(undef, 0, 11)
    skk_spillage = Array{Float64}(undef, 0, 11)
    for year = 2004:2020
        println(year)
        mat = readmatrix("$folder/Vattenf_skeälven_$year.xlsx", "Blad1", "X", year)
        skk_discharge = vcat(skk_discharge, mat[:, [1;4;6;8;10;12;14;16;18;20;22]])
        skk_spillage = vcat(skk_spillage, mat[:, [2;5;7;9;11;13;15;17;19;21;23]])
    end

    println("Reading Vattenfall production data...")
    file1 = "Skellefteälven Bruttovärden (normaltid)  2004-2018 - Vattenfall.xlsx"
    file2 = "Skellefteälven Bruttovärden (normaltid) 2019 - Vattenfall.xlsx"
    vf_raw1 = XLSX.readdata("$folder/$file1", "Rubriker", "C1:L131499")
    vf_raw2 = XLSX.readdata("$folder/$file2", "Rubriker", "C1:L8763")
    #vf_headers = string.(vf_raw1[1:3,:])
    vf = float.([vf_raw1[4:end,:]; vf_raw2[4:end,:]])                   # 2004:2019
    vattenfall_prod = [vf[:,1]  vf[:,3]+vf[:,5]  vf[:,7]+vf[:,9]]       # Bastusel, Gallejaur, Vargfors
    vattenfall_discharge = [vf[:,2]  vf[:,4]+vf[:,6]  vf[:,8]+vf[:,10]] # Bastusel, Gallejaur, Vargfors

    println("Reading Vargfors reservoir levels...")
    file = "Vargfors_vattenytor (normaltid)  2004.01.01-2019.12.31.xlsx"
    levels_vargfors = readmatrix("$folder/$file", "Rubriker", "E4:F140259")[:,[2,1]]     # 2000:2019

    println("Reading Bastusel reservoir levels...")
    file = "Bastusel_vattenytor (normaltid)  2004.01.01-2019.12.31.xlsx"
    levels_bastusel = readmatrix("$folder/$file", "Rubriker", "D4:E140259")[:,[2,1]]     # 2000:2019

    # Notera att för Gallejaure skickar vi två ytor, Gallejaure och Sandfors för säkerhets skull.
    # Gallejaure är ett lite speciellt magasin som egentligen består av två delar som vid höga vattenstånd
    # kan bli ett men vid lägre är det fallförluster mellan dessa. Kanske detta är mer än vad dina beräkningar
    # kräver men för säkerhets skull får du bägge.
    println("Reading Gallejaur reservoir levels...")
    file = "Gallejaur_vattenytor (normaltid)  2004.01.01-2019.12.31.xlsx"
    levels_gallejaur = readmatrix("$folder/$file", "Rubriker", "E4:G140259")[:,[3,1]]     # 2000:2019
    
    println("Reading Kvistforsen production data...")
    prod_kvist = Float64[]
    sheets = ["2004-2010", "2011-2016", "2017-2019"]
    lastrows = [61376, 52616, 26288]
    for i = 1:3
        sheet, lastrow = sheets[i], lastrows[i]
        println(sheet)
        mat = readmatrix("$folder/Kvistforsen MWh per h till Hanna - Chalmers - ÖVYNVY.xlsx", sheet, "B9:B$lastrow")
        prod_kvist = vcat(prod_kvist, mat)
    end

    println("Reading other Kvistforsen data...")
    kvist = Array{Float64}(undef, 0, 6)
    sheets = ["04-10 Q H W", "11-17 Q H W", "18-19 Q H W"]
    lastrows = [61369, 61369, 17528]
    columns = [[1,2,3,4,7,8], [1,2,3,4,5,6], [1,2,3,4,5,6]]
    for i = 1:3
        sheet, lastrow, cols = sheets[i], lastrows[i], columns[i]
        firstrow = (i == 3) ? 9 : 2
        println(sheet)
        mat = readmatrix("$folder/Kvistforsen MWh per h till Hanna - Chalmers - ÖVYNVY.xlsx",
                    sheet, "B$firstrow:I$lastrow"; swap=" "=>-Inf)
        kvist = vcat(kvist, mat[:,cols])
    end
    discharge_kvist, spill_kvist, upper_bergsbyn, upper_kvist, lower_kvist =
            kvist[:,2], kvist[:,3].-kvist[:,2], kvist[:,4], kvist[:,5], kvist[:,6]

    println("Reading daily water levels in several lakes...")    
    lakes_daily = readmatrix("$folder/../Tillrinning+Vattenstånd 2009-2018 (10 älvar)/Skellefteälven historiska data 2009-2018_v2 plusmeny.xlsx",
                    "Vattenstånd", "B2:G4384")[:, [2,1,3,5,6]]
                    # 2009:2018 (rader = dagar, kolumner = [Sädvajaure, Rebnisjaure, Hornavan, Storavan-Uddjaur, Naustajaure])
                    # untill 31/12/2021 for Hornavan
    println("Reading weekly reservoir levels...")
    reservoir_weekly = readmatrix("$folder/6_Totala Magasinet Standardvecka 1986-2020 med ny DGvolym.xlsx",
                    "Totala magsinet i DE", "B2:AJ53")  # 1986:2020 (rader = veckor, kolumner = år)

    println("Reading weekly inflow...")
    #file = "Tillrinningar Års- Månads- och Veckomedel och dygnsvärden med max och min 2000-2019 till Hanna_Ek_Falt.xlsx"
    #inflow_raw = XLSX.readdata("$folder/$file", "Grunddata veckomedeltillr", "D2:AG1042")
    #inflow_headers, inflow_weekly = string.(inflow_raw[1,:]), float.(inflow_raw[2:end,:]) # 2000:2019 (*52)
    file = "Skellefteälven HBV.xlsx"   # 1962-01-01:2020-09-30
    inflow_raw = XLSX.readdata("$folder/$file", "Lokaltillrinningar", "B13881:R21185")
    inflow_daily = float.(inflow_raw[:,[2;1;3:end]]) # 2000-01-01:2019-12-31

    #inflow = [inflow_weekly[:,1:5] diff(inflow_weekly[:,5:16], dims=2) diff(inflow_weekly[:,15:16], dims=2)]  #antar bergsbydammen samma som kvistforsen, borde inte spela nån roll pga ingen spillbegränsning
    #inflow_daily = vcat([repeat(inflow[y*52 + w, :]', outer=(w==52 || (w==9 && isleapyear(2000+y))) ? 8 : 7)
    #                        for y=0:19 for w=1:52]...)

    
    hours = DateTime(2004):Hour(1):DateTime(2021)-Hour(1)
    year = Dates.year.(hours)

    nhours = length(hours)
    nvf = size(vattenfall_prod, 1)

    production = fill(-Inf, nhours, 17)
    forebay_level = copy(production)
    tail_level = copy(production)
    discharge = copy(production)
    spillage = copy(production)
    lake_levels = fill(-Inf, nhours, 5)

    plants = riverdata(:Skellefte)[1]
    ix = Dict(p.name => i for (i,p) in enumerate(plants))

    skk = lookup(ix, [:Rebnis, :Sädva, :Bergnäs, :Slagnäs, :Grytfors, :Rengård,
            :Båtfors, :Finnfors, :Granfors, :Krångfors, :Selsfors])     # Skelleftekraft
    vf = lookup(ix, [:Bastusel, :Gallejaur, :Vargfors])                 # Vattenfall
    stk = [ix[:Kvistforsen]]                                            # Statkraft
    ghost = lookup(ix, [:Hornavan, :Bergsby])                           # Magasin utan kraftverk    

    production[:, skk] = prod
    production[1:nvf, vf] = vattenfall_prod
    production[1:nvf, stk] = prod_kvist
    production[1:nvf, ghost] .= 0

    forebay_level[:, skk] = level_upper
    forebay_level[1:nvf, ix[:Bastusel]] = levels_bastusel[:,1]
    forebay_level[1:nvf, ix[:Gallejaur]] = levels_gallejaur[:,1]
    forebay_level[1:nvf, ix[:Vargfors]] = levels_vargfors[:,1]
    forebay_level[1:nvf, ix[:Kvistforsen]] = upper_kvist
    forebay_level[1:nvf, ix[:Bergsby]] = upper_bergsbyn

    tail_level[:, skk] = level_lower
    tail_level[1:nvf, ix[:Bastusel]] = levels_bastusel[:,2]
    tail_level[1:nvf, ix[:Gallejaur]] = levels_gallejaur[:,2]
    tail_level[1:nvf, ix[:Vargfors]] = levels_vargfors[:,2]
    tail_level[1:nvf, ix[:Kvistforsen]] = lower_kvist
    tail_level[1:nvf, ix[:Bergsby]] = upper_bergsbyn

    day1 = sum(daysinyear.(2004:2008))*24
    lake_levels[day1+1:day1+sum(daysinyear.(2009:2020))*24, :] = repeat(lakes_daily, inner=(24,1))
    forebay_level[:, 3] = lake_levels[:, 3]
    tail_level[:, 3] = lake_levels[:, 3]

    discharge[:, skk] = skk_discharge
    discharge[1:nvf, vf] = vattenfall_discharge
    discharge[1:nvf, stk] = discharge_kvist
    discharge[1:nvf, ghost] .= 0

    spillage[:, skk] = skk_spillage
    spillage[1:nvf, vf] .= NaN
    spillage[1:nvf, stk] = spill_kvist
    spillage[1:nvf, ghost] .= NaN

    

    jldsave("$DATAFOLDER/skelleftedata_raw.jld2"; hours, production, forebay_level, tail_level,
        discharge, spillage, reservoir_weekly, inflow_daily, elspot_price, lake_levels, compress=true)
    nothing
    
end
