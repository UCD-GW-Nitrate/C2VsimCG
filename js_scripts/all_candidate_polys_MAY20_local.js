var all_candidate_polys = [
	{ id: 877, DivND: 2, Polygon: [[36.87241,-120.0439],[36.89088,-119.9575],[36.83604,-119.9202],[36.81764,-120.0133]]},
	{ id: 878, DivND: 2, Polygon: [[36.89088,-119.9575],[36.92452,-119.8448],[36.84594,-119.8309],[36.83604,-119.9202]]},
	{ id: 879, DivND: 2, Polygon: [[36.92452,-119.8448],[36.90196,-119.7709],[36.84594,-119.8309]]},
	{ id: 880, DivND: 2, Polygon: [[36.92452,-119.8448],[36.98942,-119.8064],[36.97604,-119.7281],[36.90196,-119.7709]]},
	{ id: 890, DivND: 4, Polygon: [[36.80482,-120.1214],[36.81764,-120.0133],[36.77766,-119.9811],[36.75828,-120.0732]]},
	{ id: 891, DivND: 2, Polygon: [[36.81764,-120.0133],[36.83604,-119.9202],[36.79866,-119.8683],[36.77766,-119.9811]]},
	{ id: 892, DivND: 2, Polygon: [[36.83604,-119.9202],[36.84594,-119.8309],[36.81656,-119.8129],[36.79866,-119.8683]]},
	{ id: 893, DivND: 2, Polygon: [[36.84594,-119.8309],[36.90196,-119.7709],[36.84049,-119.7426],[36.81656,-119.8129]]},
	{ id: 894, DivND: 2, Polygon: [[36.90196,-119.7709],[36.91084,-119.6954],[36.84049,-119.7426]]},
	{ id: 895, DivND: 2, Polygon: [[36.90196,-119.7709],[36.97604,-119.7281],[36.91084,-119.6954]]},
	{ id: 911, DivND: 4, Polygon: [[36.75828,-120.0732],[36.77766,-119.9811],[36.73013,-119.9317],[36.71322,-120.0211]]},
	{ id: 912, DivND: 4, Polygon: [[36.77766,-119.9811],[36.79866,-119.8683],[36.75176,-119.828],[36.73013,-119.9317]]},
	{ id: 913, DivND: 4, Polygon: [[36.79866,-119.8683],[36.81656,-119.8129],[36.7716,-119.7728],[36.75176,-119.828]]},
	{ id: 916, DivND: 2, Polygon: [[36.91084,-119.6954],[36.88117,-119.622],[36.84024,-119.6664]]},
	{ id: 926, DivND: 4, Polygon: [[36.71322,-120.0211],[36.73013,-119.9317],[36.65318,-119.9055],[36.65452,-119.985]]},
	{ id: 927, DivND: 4, Polygon: [[36.73013,-119.9317],[36.66986,-119.8459],[36.65318,-119.9055]]},
	{ id: 928, DivND: 4, Polygon: [[36.73013,-119.9317],[36.75176,-119.828],[36.69697,-119.7693],[36.66986,-119.8459]]},
	{ id: 930, DivND: 4, Polygon: [[36.7716,-119.7728],[36.79915,-119.7098],[36.75087,-119.6679],[36.7184,-119.7184]]},
	{ id: 931, DivND: 4, Polygon: [[36.79915,-119.7098],[36.84024,-119.6664],[36.78397,-119.6313],[36.75087,-119.6679]]},
	{ id: 932, DivND: 4, Polygon: [[36.84024,-119.6664],[36.88117,-119.622],[36.82302,-119.5828],[36.78397,-119.6313]]},
	{ id: 933, DivND: 2, Polygon: [[36.88117,-119.622],[36.8606,-119.5396],[36.82302,-119.5828]]},
	{ id: 967, DivND: 4, Polygon: [[36.66986,-119.8459],[36.69697,-119.7693],[36.64056,-119.699],[36.62022,-119.7986]]},
	{ id: 968, DivND: 4, Polygon: [[36.69697,-119.7693],[36.7184,-119.7184],[36.67763,-119.673],[36.64056,-119.699]]},
	{ id: 969, DivND: 4, Polygon: [[36.7184,-119.7184],[36.75087,-119.6679],[36.70758,-119.6275],[36.67763,-119.673]]},
	{ id: 970, DivND: 2, Polygon: [[36.75087,-119.6679],[36.78397,-119.6313],[36.717,-119.566],[36.70758,-119.6275]]},
	{ id: 971, DivND: 2, Polygon: [[36.78397,-119.6313],[36.82302,-119.5828],[36.76992,-119.5431],[36.717,-119.566]]},
	{ id: 972, DivND: 2, Polygon: [[36.82302,-119.5828],[36.8606,-119.5396],[36.80895,-119.5062],[36.76992,-119.5431]]},
	{ id: 982, DivND: 4, Polygon: [[36.60334,-119.8747],[36.62022,-119.7986],[36.57557,-119.7493],[36.55233,-119.8181]]},
	{ id: 983, DivND: 4, Polygon: [[36.62022,-119.7986],[36.64056,-119.699],[36.59766,-119.6518],[36.57557,-119.7493]]},
	{ id: 984, DivND: 4, Polygon: [[36.64056,-119.699],[36.67763,-119.673],[36.63393,-119.6054],[36.59766,-119.6518]]},
	{ id: 985, DivND: 4, Polygon: [[36.67763,-119.673],[36.70758,-119.6275],[36.63393,-119.6054]]},
	{ id: 986, DivND: 4, Polygon: [[36.70758,-119.6275],[36.66087,-119.5291],[36.59919,-119.5547],[36.63393,-119.6054]]},
	{ id: 987, DivND: 4, Polygon: [[36.70758,-119.6275],[36.717,-119.566],[36.66087,-119.5291]]},
	{ id: 988, DivND: 4, Polygon: [[36.717,-119.566],[36.73142,-119.4748],[36.66087,-119.5291]]},
	{ id: 989, DivND: 2, Polygon: [[36.717,-119.566],[36.76992,-119.5431],[36.73142,-119.4748]]},
	{ id: 990, DivND: 2, Polygon: [[36.76992,-119.5431],[36.80895,-119.5062],[36.78504,-119.4148],[36.73142,-119.4748]]},
	{ id: 999, DivND: 4, Polygon: [[36.52749,-119.8815],[36.55233,-119.8181],[36.49255,-119.7675]]},
	{ id: 1000, DivND: 4, Polygon: [[36.55233,-119.8181],[36.57557,-119.7493],[36.51517,-119.7045],[36.49255,-119.7675]]},
	{ id: 1001, DivND: 4, Polygon: [[36.57557,-119.7493],[36.59766,-119.6518],[36.55474,-119.6163],[36.51517,-119.7045]]},
	{ id: 1002, DivND: 4, Polygon: [[36.59766,-119.6518],[36.63393,-119.6054],[36.59919,-119.5547],[36.55474,-119.6163]]},
	{ id: 1014, DivND: 5, Polygon: [[36.44681,-119.893],[36.47591,-119.8334],[36.4361,-119.8044]]},
	{ id: 1015, DivND: 5, Polygon: [[36.47591,-119.8334],[36.49255,-119.7675],[36.46441,-119.7317],[36.4361,-119.8044]]},
	{ id: 1017, DivND: 4, Polygon: [[36.51517,-119.7045],[36.55474,-119.6163],[36.50859,-119.5737],[36.47808,-119.6223]]},
	{ id: 1018, DivND: 4, Polygon: [[36.55474,-119.6163],[36.59919,-119.5547],[36.51305,-119.5083],[36.50859,-119.5737]]},
	{ id: 1019, DivND: 4, Polygon: [[36.59919,-119.5547],[36.57236,-119.4681],[36.51305,-119.5083]]},
	{ id: 1020, DivND: 4, Polygon: [[36.59919,-119.5547],[36.66087,-119.5291],[36.63891,-119.471],[36.57236,-119.4681]]},
	{ id: 1021, DivND: 4, Polygon: [[36.66087,-119.5291],[36.73142,-119.4748],[36.63891,-119.471]]},
	{ id: 1022, DivND: 2, Polygon: [[36.73142,-119.4748],[36.71465,-119.4113],[36.63891,-119.471]]},
	{ id: 1023, DivND: 2, Polygon: [[36.73142,-119.4748],[36.78504,-119.4148],[36.74308,-119.3644],[36.71465,-119.4113]]},
	{ id: 1031, DivND: 5, Polygon: [[36.4361,-119.8044],[36.46441,-119.7317],[36.41521,-119.6846],[36.38162,-119.7661]]},
	{ id: 1041, DivND: 5, Polygon: [[36.34469,-119.8201],[36.38162,-119.7661],[36.31246,-119.7173],[36.29367,-119.7845]]},
	{ id: 1042, DivND: 5, Polygon: [[36.38162,-119.7661],[36.41521,-119.6846],[36.35759,-119.6303],[36.31246,-119.7173]]},
	{ id: 1043, DivND: 5, Polygon: [[36.41521,-119.6846],[36.42167,-119.6074],[36.38989,-119.5473],[36.35759,-119.6303]]},
	{ id: 1044, DivND: 5, Polygon: [[36.42167,-119.6074],[36.44097,-119.5572],[36.38989,-119.5473]]},
	{ id: 1045, DivND: 4, Polygon: [[36.44097,-119.5572],[36.51305,-119.5083],[36.49034,-119.4414],[36.40127,-119.4876]]},
	{ id: 1046, DivND: 4, Polygon: [[36.51305,-119.5083],[36.57236,-119.4681],[36.55996,-119.408],[36.49034,-119.4414]]},
	{ id: 1047, DivND: 4, Polygon: [[36.57236,-119.4681],[36.63891,-119.471],[36.61443,-119.3996],[36.55996,-119.408]]},
	{ id: 1048, DivND: 4, Polygon: [[36.63891,-119.471],[36.71465,-119.4113],[36.66204,-119.3853],[36.61443,-119.3996]]},
	{ id: 1049, DivND: 2, Polygon: [[36.71465,-119.4113],[36.74308,-119.3644],[36.70594,-119.3405],[36.66204,-119.3853]]},
	{ id: 1056, DivND: 5, Polygon: [[36.29367,-119.7845],[36.31246,-119.7173],[36.25992,-119.679],[36.23628,-119.7582]]},
	{ id: 1057, DivND: 5, Polygon: [[36.31246,-119.7173],[36.35759,-119.6303],[36.27986,-119.5946],[36.25992,-119.679]]},
	{ id: 1058, DivND: 5, Polygon: [[36.35759,-119.6303],[36.38989,-119.5473],[36.32612,-119.5485],[36.27986,-119.5946]]},
	{ id: 1061, DivND: 4, Polygon: [[36.40127,-119.4876],[36.49034,-119.4414],[36.45788,-119.3963],[36.41445,-119.4113]]},
	{ id: 1062, DivND: 4, Polygon: [[36.49034,-119.4414],[36.55996,-119.408],[36.51254,-119.3561],[36.45788,-119.3963]]},
	{ id: 1063, DivND: 4, Polygon: [[36.55996,-119.408],[36.61443,-119.3996],[36.58127,-119.3402],[36.51254,-119.3561]]},
	{ id: 1064, DivND: 2, Polygon: [[36.61443,-119.3996],[36.66204,-119.3853],[36.65933,-119.3007],[36.58127,-119.3402]]},
	{ id: 1065, DivND: 2, Polygon: [[36.66204,-119.3853],[36.70594,-119.3405],[36.65933,-119.3007]]},
	{ id: 1076, DivND: 6, Polygon: [[36.40127,-119.4876],[36.41445,-119.4113],[36.34858,-119.3783],[36.33617,-119.4825]]},
	{ id: 1077, DivND: 4, Polygon: [[36.51254,-119.3561],[36.4849,-119.2966],[36.45788,-119.3963]]},
	{ id: 1079, DivND: 4, Polygon: [[36.51254,-119.3561],[36.58127,-119.3402],[36.55334,-119.2716],[36.4849,-119.2966]]},
	{ id: 1080, DivND: 2, Polygon: [[36.65933,-119.3007],[36.60586,-119.2457],[36.55334,-119.2716],[36.58127,-119.3402]]},
	{ id: 1087, DivND: 6, Polygon: [[36.41445,-119.4113],[36.42238,-119.343],[36.35008,-119.3068],[36.34858,-119.3783]]},
	{ id: 1088, DivND: 6, Polygon: [[36.42238,-119.343],[36.41871,-119.2767],[36.35786,-119.2431],[36.35008,-119.3068]]},
	{ id: 1089, DivND: 4, Polygon: [[36.42238,-119.343],[36.4849,-119.2966],[36.41871,-119.2767]]},
	{ id: 1090, DivND: 4, Polygon: [[36.4849,-119.2966],[36.55334,-119.2716],[36.4839,-119.2175],[36.41871,-119.2767]]},
	{ id: 1091, DivND: 2, Polygon: [[36.55334,-119.2716],[36.60586,-119.2457],[36.54379,-119.1915],[36.4839,-119.2175]]},
	{ id: 1093, DivND: 2, Polygon: [[36.4839,-119.2175],[36.54379,-119.1915],[36.45724,-119.1091],[36.43731,-119.1698]]},
	{ id: 1094, DivND: 6, Polygon: [[36.34858,-119.3783],[36.35008,-119.3068],[36.30022,-119.2853],[36.28846,-119.3634]]},
	{ id: 1095, DivND: 6, Polygon: [[36.35008,-119.3068],[36.35786,-119.2431],[36.30776,-119.2],[36.30022,-119.2853]]},
	{ id: 1096, DivND: 6, Polygon: [[36.41871,-119.2767],[36.39551,-119.1668],[36.35694,-119.1213],[36.35786,-119.2431]]},
	{ id: 1097, DivND: 2, Polygon: [[36.43731,-119.1698],[36.45724,-119.1091],[36.4251,-119.0712],[36.39551,-119.1668]]},
	{ id: 1099, DivND: 2, Polygon: [[36.39551,-119.1668],[36.4251,-119.0712],[36.39351,-119.0327],[36.35694,-119.1213]]},
	{ id: 1101, DivND: 6, Polygon: [[36.28846,-119.3634],[36.30022,-119.2853],[36.25446,-119.2698],[36.23636,-119.3659]]},
	{ id: 1102, DivND: 6, Polygon: [[36.30022,-119.2853],[36.30776,-119.2],[36.25773,-119.1983],[36.25446,-119.2698]]},
	{ id: 1103, DivND: 6, Polygon: [[36.30776,-119.2],[36.35694,-119.1213],[36.26365,-119.1328],[36.25773,-119.1983]]},
	{ id: 1104, DivND: 2, Polygon: [[36.35694,-119.1213],[36.32273,-119.0588],[36.27668,-119.07],[36.26365,-119.1328]]},
	{ id: 1105, DivND: 2, Polygon: [[36.35694,-119.1213],[36.39351,-119.0327],[36.32273,-119.0588]]},
	{ id: 1112, DivND: 7, Polygon: [[36.1431,-119.633],[36.13709,-119.5392],[36.08321,-119.5408],[36.0865,-119.6345]]},
	{ id: 1113, DivND: 7, Polygon: [[36.1431,-119.633],[36.19114,-119.537],[36.13709,-119.5392]]},
	{ id: 1116, DivND: 6, Polygon: [[36.23636,-119.3659],[36.25446,-119.2698],[36.20795,-119.274],[36.18627,-119.3583]]},
	{ id: 1117, DivND: 6, Polygon: [[36.25446,-119.2698],[36.25773,-119.1983],[36.16689,-119.1819],[36.20795,-119.274]]},
	{ id: 1118, DivND: 6, Polygon: [[36.25773,-119.1983],[36.26365,-119.1328],[36.19979,-119.1348],[36.16689,-119.1819]]},
	{ id: 1119, DivND: 2, Polygon: [[36.26365,-119.1328],[36.27668,-119.07],[36.1976,-119.0891],[36.19979,-119.1348]]},
	{ id: 1120, DivND: 2, Polygon: [[36.27668,-119.07],[36.20604,-119.0213],[36.1976,-119.0891]]},
	{ id: 1121, DivND: 6, Polygon: [[36.13709,-119.5392],[36.16029,-119.44],[36.07919,-119.4349],[36.08321,-119.5408]]},
	{ id: 1122, DivND: 6, Polygon: [[36.16029,-119.44],[36.18627,-119.3583],[36.12666,-119.3787],[36.07919,-119.4349]]},
	{ id: 1123, DivND: 6, Polygon: [[36.18627,-119.3583],[36.13951,-119.3286],[36.12666,-119.3787]]},
	{ id: 1124, DivND: 6, Polygon: [[36.18627,-119.3583],[36.20795,-119.274],[36.15335,-119.2742],[36.13951,-119.3286]]},
	{ id: 1125, DivND: 6, Polygon: [[36.20795,-119.274],[36.16689,-119.1819],[36.15335,-119.2742]]},
	{ id: 1127, DivND: 2, Polygon: [[36.1976,-119.0891],[36.20604,-119.0213],[36.11493,-119.0143],[36.11711,-119.0621]]},
	{ id: 1129, DivND: 8, Polygon: [[36.16689,-119.1819],[36.11711,-119.0621],[36.09928,-119.1574]]},
	{ id: 1130, DivND: 8, Polygon: [[36.15335,-119.2742],[36.16689,-119.1819],[36.09928,-119.1574],[36.11935,-119.2516]]},
	{ id: 1131, DivND: 8, Polygon: [[36.13951,-119.3286],[36.15335,-119.2742],[36.11935,-119.2516],[36.1051,-119.322]]},
	{ id: 1132, DivND: 8, Polygon: [[36.12666,-119.3787],[36.13951,-119.3286],[36.1051,-119.322],[36.08749,-119.3713]]},
	{ id: 1133, DivND: 8, Polygon: [[36.07919,-119.4349],[36.12666,-119.3787],[36.08749,-119.3713]]},
	{ id: 1134, DivND: 8, Polygon: [[36.08321,-119.5408],[36.07919,-119.4349],[36.03049,-119.5409]]},
	{ id: 1135, DivND: 7, Polygon: [[36.0865,-119.6345],[36.08321,-119.5408],[36.03049,-119.5409],[36.03569,-119.6375]]},
	{ id: 1146, DivND: 7, Polygon: [[36.03569,-119.6375],[36.03049,-119.5409],[35.97834,-119.5449],[35.97301,-119.6244]]},
	{ id: 1148, DivND: 8, Polygon: [[36.03049,-119.5409],[36.07919,-119.4349],[36.02478,-119.4294],[35.99307,-119.479]]},
	{ id: 1149, DivND: 8, Polygon: [[36.07919,-119.4349],[36.08749,-119.3713],[36.03478,-119.3755],[36.02478,-119.4294]]},
	{ id: 1150, DivND: 8, Polygon: [[36.08749,-119.3713],[36.1051,-119.322],[36.04142,-119.3156],[36.03478,-119.3755]]},
	{ id: 1151, DivND: 8, Polygon: [[36.1051,-119.322],[36.11935,-119.2516],[36.05659,-119.2517],[36.04142,-119.3156]]},
	{ id: 1152, DivND: 8, Polygon: [[36.11935,-119.2516],[36.09928,-119.1574],[36.05312,-119.1613],[36.05659,-119.2517]]},
	{ id: 1153, DivND: 8, Polygon: [[36.09928,-119.1574],[36.06408,-119.1028],[36.02707,-119.1158],[36.05312,-119.1613]]},
	{ id: 1154, DivND: 8, Polygon: [[36.09928,-119.1574],[36.11711,-119.0621],[36.06408,-119.1028]]},
	{ id: 1155, DivND: 2, Polygon: [[36.11711,-119.0621],[36.11493,-119.0143],[36.04602,-119.0306],[36.06408,-119.1028]]},
	{ id: 1157, DivND: 2, Polygon: [[36.06408,-119.1028],[36.04602,-119.0306],[35.99588,-119.05],[36.02707,-119.1158]]},
	{ id: 1168, DivND: 8, Polygon: [[36.03478,-119.3755],[36.04142,-119.3156],[35.9808,-119.2989],[35.98111,-119.3652]]},
	{ id: 1169, DivND: 8, Polygon: [[36.04142,-119.3156],[36.05659,-119.2517],[35.98872,-119.2402],[35.9808,-119.2989]]},
	{ id: 1170, DivND: 8, Polygon: [[36.05659,-119.2517],[36.05312,-119.1613],[35.9895,-119.1706],[35.98872,-119.2402]]},
	{ id: 1172, DivND: 2, Polygon: [[36.02707,-119.1158],[35.99588,-119.05],[35.93288,-119.0599],[35.96683,-119.108]]},
	{ id: 1177, DivND: 2, Polygon: [[35.9895,-119.1706],[35.96683,-119.108],[35.89791,-119.1483],[35.92973,-119.1887]]},
	{ id: 1183, DivND: 2, Polygon: [[35.92973,-119.1887],[35.89791,-119.1483],[35.82683,-119.159],[35.87191,-119.2175]]},
	{ id: 1203, DivND: 2, Polygon: [[35.87191,-119.2175],[35.82683,-119.159],[35.8256,-119.2451]]},
	{ id: 1208, DivND: 2, Polygon: [[35.8256,-119.2451],[35.82683,-119.159],[35.77732,-119.1549],[35.77735,-119.2637]]},
	{ id: 1220, DivND: 2, Polygon: [[35.78,-119.4688],[35.77946,-119.4051],[35.70029,-119.4068],[35.7133,-119.472]]},
	{ id: 1221, DivND: 2, Polygon: [[35.77946,-119.4051],[35.7802,-119.3543],[35.70969,-119.3379],[35.70029,-119.4068]]},
	{ id: 1223, DivND: 2, Polygon: [[35.77735,-119.2637],[35.77732,-119.1549],[35.72342,-119.1798],[35.71778,-119.266]]},
	{ id: 1232, DivND: 3, Polygon: [[35.72065,-119.6534],[35.72397,-119.5843],[35.66932,-119.5972],[35.66756,-119.6775]]},
	{ id: 1233, DivND: 3, Polygon: [[35.72397,-119.5843],[35.72501,-119.5215],[35.6662,-119.5331],[35.66932,-119.5972]]},
	{ id: 1235, DivND: 2, Polygon: [[35.7133,-119.472],[35.70029,-119.4068],[35.64648,-119.4056],[35.65043,-119.4751]]},
	{ id: 1236, DivND: 2, Polygon: [[35.70029,-119.4068],[35.70969,-119.3379],[35.65069,-119.3295],[35.64648,-119.4056]]},
	{ id: 1238, DivND: 2, Polygon: [[35.71778,-119.266],[35.72342,-119.1798],[35.66283,-119.1752],[35.65644,-119.2625]]},
	{ id: 1247, DivND: 3, Polygon: [[35.66756,-119.6775],[35.66932,-119.5972],[35.62095,-119.5969],[35.60938,-119.6753]]},
	{ id: 1248, DivND: 3, Polygon: [[35.66932,-119.5972],[35.6662,-119.5331],[35.62167,-119.534],[35.62095,-119.5969]]},
	{ id: 1249, DivND: 2, Polygon: [[35.6662,-119.5331],[35.65043,-119.4751],[35.5909,-119.4746],[35.62167,-119.534]]},
	{ id: 1251, DivND: 2, Polygon: [[35.64648,-119.4056],[35.65069,-119.3295],[35.59788,-119.4126]]},
	{ id: 1252, DivND: 2, Polygon: [[35.65069,-119.3295],[35.59758,-119.3205],[35.59788,-119.4126]]},
	{ id: 1253, DivND: 3, Polygon: [[35.65069,-119.3295],[35.65644,-119.2625],[35.5971,-119.2471],[35.59758,-119.3205]]},
	{ id: 1254, DivND: 2, Polygon: [[35.65644,-119.2625],[35.66283,-119.1752],[35.6009,-119.1689],[35.5971,-119.2471]]},
	{ id: 1265, DivND: 2, Polygon: [[35.62095,-119.5969],[35.62167,-119.534],[35.5909,-119.4746],[35.57688,-119.5739]]},
	{ id: 1270, DivND: 2, Polygon: [[35.57688,-119.5739],[35.5909,-119.4746],[35.53369,-119.4525],[35.52638,-119.5394]]},
	{ id: 1271, DivND: 2, Polygon: [[35.5909,-119.4746],[35.59788,-119.4126],[35.53387,-119.3898],[35.53369,-119.4525]]},
	{ id: 1272, DivND: 2, Polygon: [[35.59788,-119.4126],[35.59758,-119.3205],[35.53167,-119.3156],[35.53387,-119.3898]]},
	{ id: 1273, DivND: 2, Polygon: [[35.59758,-119.3205],[35.5971,-119.2471],[35.53232,-119.2409],[35.53167,-119.3156]]},
	{ id: 1274, DivND: 2, Polygon: [[35.5971,-119.2471],[35.6009,-119.1689],[35.54302,-119.159],[35.53232,-119.2409]]},
	{ id: 1284, DivND: 3, Polygon: [[35.53387,-119.3898],[35.53167,-119.3156],[35.48229,-119.3045],[35.4801,-119.3669]]},
	{ id: 1285, DivND: 3, Polygon: [[35.53167,-119.3156],[35.53232,-119.2409],[35.48261,-119.238],[35.48229,-119.3045]]},
	{ id: 1286, DivND: 2, Polygon: [[35.53232,-119.2409],[35.54302,-119.159],[35.491,-119.1634],[35.48261,-119.238]]},
	{ id: 1287, DivND: 3, Polygon: [[35.54302,-119.159],[35.5531,-119.0915],[35.49951,-119.1102],[35.491,-119.1634]]},
	{ id: 1296, DivND: 3, Polygon: [[35.4801,-119.3669],[35.44202,-119.2938],[35.42089,-119.3732]]},
	{ id: 1298, DivND: 3, Polygon: [[35.48229,-119.3045],[35.48261,-119.238],[35.43222,-119.2391],[35.44202,-119.2938]]},
	{ id: 1299, DivND: 2, Polygon: [[35.48261,-119.238],[35.491,-119.1634],[35.43168,-119.1486],[35.43222,-119.2391]]},
	{ id: 1300, DivND: 2, Polygon: [[35.491,-119.1634],[35.49951,-119.1102],[35.46428,-119.0794],[35.43168,-119.1486]]},
	{ id: 1308, DivND: 3, Polygon: [[35.42089,-119.3732],[35.44202,-119.2938],[35.38488,-119.2923],[35.38164,-119.3696]]},
	{ id: 1309, DivND: 3, Polygon: [[35.44202,-119.2938],[35.43222,-119.2391],[35.38,-119.2226],[35.38488,-119.2923]]},
	{ id: 1311, DivND: 2, Polygon: [[35.43168,-119.1486],[35.46428,-119.0794],[35.43319,-119.0518],[35.38484,-119.1346]]},
	{ id: 1320, DivND: 3, Polygon: [[35.38164,-119.3696],[35.38488,-119.2923],[35.33827,-119.2419],[35.33109,-119.3164]]},
	{ id: 1321, DivND: 3, Polygon: [[35.38488,-119.2923],[35.38,-119.2226],[35.33364,-119.167],[35.33827,-119.2419]]},
	{ id: 1322, DivND: 3, Polygon: [[35.38,-119.2226],[35.38484,-119.1346],[35.3474,-119.0887],[35.33364,-119.167]]},
	{ id: 1323, DivND: 2, Polygon: [[35.38484,-119.1346],[35.43319,-119.0518],[35.38004,-119.0359],[35.3474,-119.0887]]},
	{ id: 1324, DivND: 3, Polygon: [[35.43319,-119.0518],[35.41,-118.9798],[35.38004,-119.0359]]},
	{ id: 1325, DivND: 3, Polygon: [[35.38004,-119.0359],[35.41,-118.9798],[35.35236,-118.934],[35.34571,-119.0012]]},
	{ id: 1328, DivND: 3, Polygon: [[35.38004,-119.0359],[35.34571,-119.0012],[35.3474,-119.0887]]},
	{ id: 1330, DivND: 3, Polygon: [[35.30882,-119.3663],[35.33109,-119.3164],[35.2921,-119.2628],[35.26053,-119.3242]]},
	{ id: 1331, DivND: 3, Polygon: [[35.33109,-119.3164],[35.33827,-119.2419],[35.30915,-119.2097],[35.2921,-119.2628]]},
	{ id: 1332, DivND: 3, Polygon: [[35.33827,-119.2419],[35.33364,-119.167],[35.30915,-119.2097]]},
	{ id: 1333, DivND: 3, Polygon: [[35.33364,-119.167],[35.3474,-119.0887],[35.28467,-119.0413],[35.28984,-119.116]]},
	{ id: 1334, DivND: 3, Polygon: [[35.3474,-119.0887],[35.34571,-119.0012],[35.28775,-118.9727],[35.28467,-119.0413]]},
	{ id: 1335, DivND: 3, Polygon: [[35.34571,-119.0012],[35.35236,-118.934],[35.28917,-118.8952],[35.28775,-118.9727]]},
	{ id: 1344, DivND: 3, Polygon: [[35.2921,-119.2628],[35.30915,-119.2097],[35.24184,-119.1299],[35.23545,-119.2087]]},
	{ id: 1345, DivND: 3, Polygon: [[35.30915,-119.2097],[35.33364,-119.167],[35.28984,-119.116],[35.24184,-119.1299]]},
	{ id: 1346, DivND: 3, Polygon: [[35.28984,-119.116],[35.28467,-119.0413],[35.2327,-119.0328],[35.24184,-119.1299]]},
	{ id: 1347, DivND: 3, Polygon: [[35.28467,-119.0413],[35.28775,-118.9727],[35.22594,-118.9601],[35.2327,-119.0328]]},
	{ id: 1348, DivND: 3, Polygon: [[35.28775,-118.9727],[35.28917,-118.8952],[35.24246,-118.8828],[35.22594,-118.9601]]},
	{ id: 1354, DivND: 3, Polygon: [[35.23545,-119.2087],[35.24184,-119.1299],[35.19688,-119.1311],[35.1758,-119.2113]]},
	{ id: 1355, DivND: 3, Polygon: [[35.24184,-119.1299],[35.2327,-119.0328],[35.16711,-119.0498],[35.19688,-119.1311]]},
	{ id: 1356, DivND: 3, Polygon: [[35.2327,-119.0328],[35.22594,-118.9601],[35.17415,-118.9734],[35.16711,-119.0498]]},
	{ id: 1357, DivND: 3, Polygon: [[35.22594,-118.9601],[35.24246,-118.8828],[35.19095,-118.8615],[35.17415,-118.9734]]},
	{ id: 1358, DivND: 2, Polygon: [[35.24246,-118.8828],[35.26896,-118.8183],[35.22901,-118.7839],[35.19095,-118.8615]]},
	{ id: 1365, DivND: 3, Polygon: [[35.1758,-119.2113],[35.19688,-119.1311],[35.12787,-119.14],[35.12783,-119.2122]]},
	{ id: 1366, DivND: 3, Polygon: [[35.19688,-119.1311],[35.16711,-119.0498],[35.08978,-119.0604],[35.12787,-119.14]]},
	{ id: 1368, DivND: 3, Polygon: [[35.17415,-118.9734],[35.19095,-118.8615],[35.13191,-118.9069],[35.10823,-118.9877]]},
	{ id: 1369, DivND: 2, Polygon: [[35.19095,-118.8615],[35.14374,-118.8255],[35.09741,-118.8884],[35.13191,-118.9069]]},
	{ id: 1370, DivND: 2, Polygon: [[35.19095,-118.8615],[35.22901,-118.7839],[35.18434,-118.7514],[35.14374,-118.8255]]},
];