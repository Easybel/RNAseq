#Here, we first paste all the 
IDwt="SS_root_wt_"
IDCu="SS_root_wt_"
end="_raw.count"
IDout="myTab"
paste -d' ' $IDwt"s1"$end $IDwt"s2"$end $IDwt"s3"$end $IDCu"s1"$end $IDCu"s2"$end $IDCu"s3"$end| grep -v "^#" | awk '{print $1,$7,$14,$21,$28,$35,$42}' | tail -n +2 > $IDout".csv"
 
