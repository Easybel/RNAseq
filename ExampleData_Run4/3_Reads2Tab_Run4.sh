#Here, we first paste all the files next to each other, then we print (awk) the columns we want 
dataPath="/home/isabel/sciebo/P5_ExpEvol_Ngo/RNA/2022March_Run4/1a_Result_MapCount/Map2_MS11Ref"

cd $dataPath

IDctl="Run4_NgoG4"
#ID1="Run4_NgoMS11"
#ID2="Run4_NgoG4pilE"
#ID3="Run4_NgoG4pilT"
ID1="Run4_NgoG4pilENG17"
ID2="Run4_NgoG4pilENG24"
ID3="Run4_NgoG4pilENG32"

end="_O_minO20_CDS_raw.count"
IDout="Run4Ngo_CDS_Counts_G4_NG17_NG24_NG32"
#paste $IDctl"_S1"$end $IDctl"_S2"$end $IDctl"_S3"$end $ID1"_S1"$end $ID1"_S2"$end $ID1"_S3"$end $ID2"_S1"$end $ID2"_S2"$end $ID2"_S3"$end $ID3"_S1"$end $ID3"_S2"$end $ID3"_S3"$end $ID4"_S1"$end $ID4"_S2"$end $ID4"_S3"$end $ID5"_S1"$end $ID5"_S2"$end $ID5"_S3"$end $ID6"_S1"$end $ID6"_S2"$end $ID6"_S3"$end | grep -v "^#" | awk '{print $1,$7,$14,$21,$28,$35,$42,$49,$56,$63,$70,$77,$84,$91,$98,$105,$112,$119,$126,$133,$140,$147}' > $IDout".csv"

paste $IDctl"_S1"$end $IDctl"_S2"$end $IDctl"_S3"$end $ID1"_S1"$end $ID1"_S2"$end $ID1"_S3"$end $ID2"_S1"$end $ID2"_S2"$end $ID2"_S3"$end $ID3"_S1"$end $ID3"_S2"$end $ID3"_S3"$end | grep -v "^#" | awk '{print $1,$7,$14,$21,$28,$35,$42,$49,$56,$63,$70,$77,$84,$91,$98,$105,$112,$119,$126,$133,$140,$147}' > $IDout".csv"


#optionally: we remove the first line (tail) 
#paste -d' ' $IDwt"s1"$end $IDwt"s2"$end $IDwt"s3"$end $IDCu"s1"$end $IDCu"s2"$end $IDCu"s3"$end| grep -v "^#" | awk '{print $1,$7,$14,$21,$28,$35,$42}' | tail -n +2 > $IDout".csv"


