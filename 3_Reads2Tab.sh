#Here, we first paste all the files next to each other, then we print (awk) the columns we want 
IDctl="Ngo_T1ctl"
IDA="Ngo_T1Azi"
end="_O_minO20_CDS_raw.count"
IDout="Ngo_CDS_Counts"
paste -d' ' $IDctl"_S1"$end $IDctl"_S2"$end $IDctl"_S3"$end $IDA"_S1"$end $IDA"_S2"$end $IDA"_S3"$end| grep -v "^#" \
| awk '{print $1,$7,$14,$21,$28,$35,$42}' > $IDout".csv"

#optionally: we remove the first line (tail) 
#paste -d' ' $IDwt"s1"$end $IDwt"s2"$end $IDwt"s3"$end $IDCu"s1"$end $IDCu"s2"$end $IDCu"s3"$end| grep -v "^#" | awk '{print $1,$7,$14,$21,$28,$35,$42}' | tail -n +2 > $IDout".csv"


