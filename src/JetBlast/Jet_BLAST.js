	////////////////////////////////////////////////////////////////// 
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//:::::::::::::::::::::: -JETBlastReport- :::::::::::::::::::::://
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//////////////////////////////////////////////////////////////////

  // -----------------------------------------------------------------------------------
  // JETBlastReport is free software; you can redistribute it and/or
  // modify it under the terms of the GNU General Public License
  // as published by the Free Software Foundation; either version 3
  // of the License, or (at your option) any later version.
  //
  // This program is distributed in the hope that it will be useful,
  // but WITHOUT ANY WARRANTY; without even the implied warranty of
  // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  // GNU General Public License for more details.
  //
  // You should have received a copy of the GNU General Public License
  // along with this program; if not, write to the Free Software
  // Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
  // or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
  // -----------------------------------------------------------------------------------
  // This file is part of the JETBlastReport local BLAST translating program,
  // written by
  //    Benoit Piegu    (CNRS, FR)                   <benoit.piegu@tours.inra.fr>
  //	  Valentin Marcon (Université d'Auvergne, FR)  <valentin.marcon@etu.udamail.fr>
  // -----------------------------------------------------------------------------------


	//////////////////////////////////////////////////////////////////	
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//::::::::::::::::::: -General Script- ::::::::::::::::::::::::://
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//////////////////////////////////////////////////////////////////



//	Exclude the utilisation of "Internet Explorer" 
//	---------------------------------------------------------------
if (navigator.appName == "Microsoft Internet Explorer" || navigator.appName == "Windows Internet Explorer") {
 	document.body.innerHTML = '<img src="http://pocketpalsson.com/wp-content/uploads/2012/07/9148d9a92**5ab5aee46d9ebdaa0eebc44.jpg"><br><br>You can\'t use the script on IE... Please <a href="http://www.mozilla.org/fr/firefox/new/">change</a> your browser.';
}
//	---------------------------------------------------------------		


else{
	document.head.innerHTML = '<link rel="shortcut icon" type="image/png" href="http://www.aht.li/2075540/icone.png" />';
	var tab_ligne = new Array();
	var report = document.documentElement.innerHTML;
	report = report.replace(/<[^>]*>/g, "");		
	report = report.replace(/&gt;/g, ">");
	tab_ligne = report.split('\n');
	var blast_hash = new Object();
	report=null;

//	Parse the Blast Report
//-------------------------------------
	blast_hash = Read_Blast(tab_ligne);

	tab_ligne=null;
//-------------------------------------


//	Htmlize the blast with the hash of result
//---------------------------------------------
	try{
		t_align=60;
		nbhsp=10000;
		scorehsp=0;
		//Default parameters for the alignement and the graph
		Do_Html();//THE function
		Hide("Visual_Parameters"); 
        document.getElementById("help_capture").style.display = "none"; 
		//For default hide the parameters
	}

	catch(e){
		document.body.innerHTML = '<img src="http://www.aht.li/2062911/JETBlastReportError.png" alt="Error"><br><br><ul><li>'+e+'</li><br><li>Dont click 2 time on the BookLet</li><li>Maybe your report cant be read by JETBlastReport</li><br><br><li>Contact <b>benoit.piegu&nbsp;@&nbsp;tours.inra.fr</b> or <b>valentin.marcon&nbsp;@&nbsp;etu.udamail.fr</b> for more information.</li></ul>';
	}
//---------------------------------------------

}

	//////////////////////////////////////////////////////////////////
 	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
 	//::::::::::::::::: -Functions to parse- ::::::::::::::::::::::://
 	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//////////////////////////////////////////////////////////////////

 function Read_Blast(tab_l) {

 	reporttype = "unknow";
 	Squery = 1;
 	Ssbjct = 1;
 	strand = 0;
 	var ligne;
 	var bool = 0; //booleens
 	var bool2 = 0;
 	var gapped_stats = 0;
 	var resultat = new Object();
	var THE_taille=0;//bug of string parse
 	resultat["Parameters_allowgaps"] = 'yes';

 	for (var i = 0; i < tab_l.length; i++) {
 		ligne = tab_l[i];

//  Test of RegExp on every line to insert the maximum of information in a hash
//------------------------------------------------------------------------------
 		if (ligne.match(/No hits found/i)) {
 			break;
 		} else if (ligne.match(/^(\s*)?$/)) {
 			continue; //Lignes vides
 		} else if (ligne.match(/^((?:\S+?)?BLAST[NPX]?)\s*(.*)/i) || ligne.match(/^(P?GENEWISE|HFRAME|SWN|TSWN)\s+(.+)/i)) {
 			reporttype = RegExp.$1;
 			resultat["BlastOutput_version"] = RegExp.$2;
 			if (reporttype == "RPS-BLAST") {
 				reporttype = reporttype + "(BLASTP)";
 			}
 		} else if (ligne.match(/^Reference:\s+(.*)$/)) {
 			var algorithm_reference = RegExp.$1;
 			var j = 1 + i;
 			var ligne2 = tab_l[j];
 			while (!ligne2.match(/^$/) && !ligne2.match(/^RID:/) && !ligne2.match(/^Database:/) && !ligne2.match(/^Query=/)) {
 				j = j + 1;
 				algorithm_reference = algorithm_reference + ligne2;
 				ligne2 = tab_l[j];
 			}
 			resultat["BlastOutput_algorithm-reference"] = algorithm_reference;
 		} else if (ligne.match(/^RID:\s+(.*)$/)) {
 			resultat["BlastOutput_rid"] = RegExp.$1;
 		} else if (ligne.match(/^Query=\s*(\S*)\s*(.*)$/)) {
 			if (!RegExp.$1.match(/^\s*$/)) {
 				resultat["BlastOutput_query-def"] = RegExp.$1;
 			}
 		} else if (ligne.match(/\((\-?[\d,]+)\s+letters.*\)/ || /^Length=(\-?[\d,]+)/)) {
 			var intermed = RegExp.$1;
 			intermed = intermed.replace(/\,/g, '');
 			resultat["BlastOutput_query-len"] = intermed;
 		} else if (ligne.match(/Length of query:\s*(\d*)/)) {
 			resultat["BlastOutput_query-len"] = RegExp.$1;
 		} else if (ligne.match(/Length=(\d*)/)) {
 			resultat["BlastOutput_query-len"] = RegExp.$1; 
			//special blast (ex: bug2246)
 		} else if (ligne.match(/^>Unfiltered[+-]1$/)) {
 			// skip all of the lines of unfiltered sequence
 			var j = i;
 			var ligne2 = tab_l[j];
 			while (!ligne2.match(/\s/)) {
 				tab_l[j] = "";
 				j = j + 1;
 				ligne2 = tab_l[j];
 			}
 		} else if (ligne.match(/Number of Sequences:\s+([\d\,]+)/i || /of sequences in database:\s+(\-?[\d,]+)/i)) {
 			resultat["Statistics_db-len"] = RegExp.$1;
 		} else if (ligne.match(/letters in database:\s+(\-?[\d,]+)/i)) {
 			resultat["Statistics_db-let"] = RegExp.$1;
 		} else if (ligne.match(/^Gapped/)) {
 			gapped_stats = 1;
 		} else if (ligne.match(/^Lambda\s+/)) {
 			bool = 1;
 		} else if (ligne.match(/(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+)/) && bool == 1 && gapped_stats == 0) {
 			resultat["Statistics_lambda"] = RegExp.$1;
 			resultat["Statistics_kappa"] = RegExp.$2;
 			resultat["Statistics_entropy"] = RegExp.$3;
 			bool = 0;
 		} else if (ligne.match(/(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+)/) && bool == 1 && gapped_stats == 1) {
 			resultat["Statistics_gapped_lambda"] = RegExp.$1;
 			resultat["Statistics_gapped_kappa"] = RegExp.$2;
 			resultat["Statistics_gapped_entropy"] = RegExp.$3;
 			bool = 0;
 			gapped_stats = 0;
 		} else if (ligne.match(/^\s+(\-?[\d,]+|\S+)\s+sequences\;\s+(\-?[\d,]+|\S+)\s+total\s+letters/)) {
 			var intermed = RegExp.$1;
 			var intermed2 = RegExp.$2;
 			intermed = intermed.replace(/\,/g, '');
 			resultat["BlastOutput_db-len"] = intermed;
 			intermed2 = intermed2.replace(/\,/g, '');
 			resultat["BlastOutput_db-let"] = intermed2;
 		}
 		else if (ligne.match(/E=(\S+)/)) {
 			resultat["Parameters_expect"] = RegExp.$1;
 		} else if (ligne.match(/nogaps/i)) {
 			resultat["Parameters_allowgaps"] = 'no';
 		} else if (ligne.match(/^\s*Histogram/)) {
 			resultat["Parameters_allowgaps"] = 'yes';
 		} else if (ligne.match(/ctxfactor=(\S+)/)) {
 			resultat["Statistics_ctxfactor"] = RegExp.$1;
 		} else if (ligne.match(/(postsw|links|span[12]?|warnings|notes|gi|noseqs|qres|qype)/)) {
 			//			resultat["Parameters_"+RegExp.$1]='yes';
 		} else if (ligne.match(/(Frame|Strand)\s+MatID\s+Matrix name/i)) {
 			var firstgapinfo = 1;
 			var frame = 0;
 			var k = i;
 			var ligne2;
 			for (var a = 0; a < 12; a++) {
 				k = k + 1;
 				if (ligne.match(/^(\s*)?$/)) {
 					break;
 				}
 				ligne2 = tab_l[k];
 				ligne2 = ligne2.replace(/^\s+/, "");
 				ligne2 = ligne2.replace(/\s+$/, "");
 				if (ligne2.match(/Q=(\d+),R=(\d+)\s+/) && firstgapinfo == 1) {
 					firstgapinfo = 0;
 					resultat["Parameters_gap-open"] = RegExp.$1;
 					resultat["Parameters_gap-extend"] = RegExp.$2;
 					var fields = ligne2.split(/\s+/);
 					resultat["Statistics_frame" + frame + "_lambda_gapped"] = fields[1];
 					resultat["Statistics_frame" + frame + "_kappa_gapped"] = fields[2];
 					resultat["Statistics_frame" + frame + "_entropy_gapped"] = fields[3];
 				} else {
 					firstgapinfo = 1;
 					var fields = ligne2.split(/\s+/);
 					if (frame == 0) {
 						resultat["Parameters_matrix"] = fields[2];
 						resultat["Statistics_lambda"] = fields[3];
 						resultat["Statistics_kappa"] = fields[4];
 						resultat["Statistics_entropy"] = fields[5];
 					}
 					frame = fields[0];
 					var ii = 3;
 					var type = "lambda_used kappa_used entropy_used lambda_computed kappa_computed entropy_computed".split(/\s+/);
 					for (var b = 0; b < type.length; b++) {
 						var f = fields[ii];
 						if (f == "same") {
 							f = fields[ii - 3];
 						}
 						ii = ii + 1;
 						resultat["Statistics_frame" + frame + "_" + type[b]] = f;
 					}
 				}
 			}
 		} else if (ligne.match(/(Frame|Strand)\s+MatID\s+Length/i)) {
 			var frame;
 			var k = i;
 			var ligne2;
 			for (var a = 0; a < 12; a++) {
 				k = k + 1;
 				ligne2 = tab_l[k];
 				if (ligne.match(/^(\s*)?$/)) {
 					break;
 				}
 				ligne2 = ligne2.replace(/^\s+/, "");
 				ligne2 = ligne2.replace(/\s+$/, "");
 				var fields = ligne2.split(/\s+/);
 				if (fields.length <= 3) {
 					var type = "X_gapped E2_gapped S2_gapped".split(/\s+/);
 					for (var b = 0; b < type.length; b++) {
 						resultat["Statistics_frame" + frame + "_" + type[b]] = fields[b];
 					}
 				} else {
 					var type = "length efflength E S W T X E2 S2".split(/\s+/);
 					for (var b = 0; b < type.length; b++) {
 						frame = fields[0];
 						resultat["Statistics_frame" + frame + "_" + type[b]] = fields[b + 2];
 					}
 				}
 			}
 		} else if (ligne.match(/(\S+)=(\S+)/)) {
 			resultat["Parameters_" + RegExp.$1] = RegExp.$2;
 		} else if (ligne.match(/(\S+\s+\S+)\s+DFA:\s+(\S+)\s+\((.+)\)/)) {
 			if (RegExp.$1 == "state in") {
 				resultat["Statistics_DFA_states"] = RegExp.$2 + RegExp.$3;
 			} else if (RegExp.$1 == "size of") {
 				resultat["Statistics_DFA_size"] = RegExp.$2 + RegExp.$3;
 			}
 		} else if (ligne.match(/^\s+Time to generate neighborhood:\s+(\S+\s+\S+\s+\S+)/)) {
 			resultat["Statistics_neighbortime"] = RegExp.$1;
 		} else if (ligne.match(/processors\s+used:\s+(\d+)/)) {
 			resultat["Statistics_noprocessors"] = RegExp.$1;
 		} else if (ligne.match(/^\s+(\S+)\s+cpu\s+time:\s+(\S+\s+\S+\s+\S+)\s+Elapsed:\s+(\S+)/)) {
 			var cputype = RegExp.$1;
 			cputype = cputype.toLowerCase();
 			resultat["Statistics_" + cputype + "_cputime"] = RegExp.$2;
 			resultat["Statistics_" + cputype + "_actualtime"] = RegExp.$3;
 		} else if (ligne.match(/^\s+Start:/)) {
 			var elem = ligne.split(/\s+(Start|End)\:\s+/);
 			resultat["Statistics_starttime"] = Chomp(elem[2]);
 			resultat["Statistics_endtime"] = elem[4];
 		} else if (ligne.match(/^\s*Database:\s*(.*)/) && bool2 != 1) {
 			bool2 = 1;
 			resultat["Parameters_full_dbpath"] = RegExp.$1;
 			resultat["Parameters_full_dbpath"] = Chomp(resultat["Parameters_full_dbpath"]);
 			ligne2 = tab_l[i + 1];
 			if (!ligne2.match(/^\s+(\-?[\d,]+|\S+)\s+sequences\;\s+(\-?[\d,]+|\S+)\s+total\s+letters/)) {
 				ligne2 = ligne2.substr(0, ligne2.length);
 				ligne2 = Chomp(ligne2);
 				resultat["Parameters_full_dbpath"] += ligne2;
 				resultat["Parameters_full_dbpath"] = Chomp(resultat["Parameters_full_dbpath"]);
 			}
 		} else if (ligne.match(/^\s+Posted\s+(.+)/)) {
 			var d = RegExp.$1;
 			resultat["Statistics_posted_date"] = Chomp(d);
 		}
 		else if (ligne.match(/^Matrix:\s+(.+)\s*$/i)) {
 			resultat["Parameters_matrix"] = RegExp.$1;
 		} else if (ligne.match(/^effective\s+search\s+space\s+used:\s+(\d+)/i)) {
 			resultat["Statistics_eff-spaceused"] = RegExp.$1;
 		} else if (ligne.match(/^effective\s+search\s+space:\s+(\d+)/i)) {
 			resultat["Statistics_eff-space"] = RegExp.$1;
 		} else if (ligne.match(/Gap\s+Penalties:\s+Existence:\s+(\d+)\,\s+Extension:\s+(\d+)/i)) {
 			resultat["Parameters_gap-open"] = RegExp.$1;
 			resultat["Parameters_gap-extend"] = RegExp.$2;
 		} else if (ligne.match(/effective\s+HSP\s+length:\s+(\d+)/)) {
 			resultat["Statistics_hsp-len"] = RegExp.$1;
 		} else if (ligne.match(/Number\s+of\s+HSP's\s+better\s+than\s+(\S+)\s+without\s+gapping:\s+(\d+)/)) {
 			resultat["Statistics_number_of_hsps_better_than_expect_value_cutoff_without_gapping"] = RegExp.$2;
 		} else if (ligne.match(/Number\s+of\s+HSP's\s+gapped:\s+(\d+)/)) {
 			resultat["Statistics_number_of_hsps_gapped"] = RegExp.$1;
 		} else if (ligne.match(/Number\s+of\s+HSP's\s+successfully\s+gapped:\s+(\d+)/)) {
 			resultat["Statistics_number_of_hsps_successfully_gapped"] = RegExp.$1;
 		} else if (ligne.match(/Length\s+adjustment:\s+(\d+)/)) {
 			resultat["Statistics_length_adjustment"] = RegExp.$1;
 		} else if (ligne.match(/(effective)\s(length)\s(of)\s(query:)\s(\d+)/i)) {
 			resultat["Statistics_query-len"] = RegExp.$5;
 			var intermed = RegExp.$1;
 			intermed = intermed.replace(/\,/g, '');
 			resultat["Statistics_query-len"] = intermed;
 		} else if (ligne.match(/(effective)\s(length)\s(of)\s(database:)\s((\d*.?\d+)+)/)) {
 			resultat["Statistics_eff-dblen"] = RegExp.$5;
 		} else if (ligne.match(/^(T|A|X1|X2|X3|S1|S2):\s+(\d+(\.\d+)?)\s+(?:\(\s*(\d+\.\d+) bits\))?/)) {
 			var v = RegExp.$2;
 			v = Chomp(v);
 			resultat["Statistics_" + RegExp.$1] = v;
 			if (RegExp.$4) {
 				resultat["Statistics_" + RegExp.$1 + "_bits"] = RegExp.$4;
 			}
 		} else if (ligne.match(/frameshift\s+window\,\s+decay\s+const:\s+(\d+)\,\s+([\.\d]+)/)) {
 			resultat["Statistics_framewindow"] = RegExp.$1;
 			resultat["Statistics_decay"] = RegExp.$2;
 		} else if (ligne.match(/^Number\s+of\s+Hits\s+to\s+DB:\s+(\S+)/)) {
 			resultat["Statistics_hit_to_db"] = RegExp.$1;
 		} else if (ligne.match(/^Number\s+of\s+extensions:\s+(\S+)/)) {
 			resultat["Statistics_num_extensions"] = RegExp.$1;
 		} else if (ligne.match(/^Number\s+of\s+successful\s+extensions:\s+(\S+)/)) {
 			resultat["Statistics_num_suc_extensions"] = RegExp.$1;
 		} else if (ligne.match(/^Number\s+of\s+sequences\s+better\s+than\s+(\S+):\s+(\d+)/)) {
 			resultat["Parameters_expect"] = RegExp.$1;
 			resultat["Statistics_seqs_better_than_cutoff"] = RegExp.$2;
 		} else if (ligne.match(/^\s+Posted\s+date:\s+(.+)/)) {
 			var d = RegExp.$1;
 			d = Chomp(d);
 			resultat["Statistics_posted_date"] = d;
 		}
 	}
 	resultat["BlastOutput_tabHit"] = Fill_Tab_hit(tab_l); //Remplir les Hits & Hsps
 	resultat["BlastOutput_program"] = reporttype; //Définir le type de blast
 	return resultat;
 }
//*************************************************************************************************************************

 function Fill_Tab_hit(tab_l) {
 	var ligne;
 	var iteration = 0;
 	var t_hit = new Object();

 	for (var i = 0; i < tab_l.length; i++) {
 		ligne = tab_l[i];
 		if (ligne.match(/^>\s*(\S+)/)) {
 			t_hit[iteration] = Extract_Hit(tab_l, i);
 			iteration = iteration + 1;
 		} else if (ligne.match(/^(\s*)?$/)) {
 			continue;
 		}
 	}
 	return t_hit;
 }
//*************************************************************************************************************************

 function Extract_Hit(tab, index) {
 	var hit = new Object();
 	var ligne;
 	var numhit = 0;
 	var i = index;
 	var gi;
 	var acc;
 	var version;
 	var bool = 0;
 	var restofline;

 	for (i; i < tab.length; i++) {
 		ligne = tab[i];
 		if (ligne.match(/^>\s*(\S+)/)) {
 			numhit = numhit + 1;
 		}
 		if (numhit <= 1) {
 			if (ligne.match(/^>\s*(\S+)\s*(.*)?/)) {
 				ind = i;
 				hit["Hit_id"] = RegExp.$1;
 				restofline = RegExp.$2;
 				acc = Get_Seq_Identifiers(RegExp.$1, ligne);
 				hit["Hit_accession"] = acc;
 				hit["Hit_def"] = restofline;
 				bool = 1;
 			} else if (ligne.match(/Length\s*=\s*([\d,]+)/)) {
 				hit["Hit_len"] = RegExp.$1;
 				bool = 0;
 			} else if (ligne.match(/^\s/) && bool == 1 && !ligne.match(/Length\s*=\s*([\d,]+)/)) {
 				restofline = restofline + ligne;
 				hit["Hit_def"] = restofline;
 			} else if (ligne.match(/^(\s*)?$/)) {
 				continue;
 			}
 		}
 		if (ligne.match(/\sScore\s=\s*([0-9]+\.?[0-9]*)/i)) {
 			break;
 		}
 		if (numhit > 1) {
 			break;
 		}
 	}
 	hit["tab_hsp"] = Fill_Tab_Hsp(tab, ind);
 	return hit;
 }
//*************************************************************************************************************************

 function Fill_Tab_Hsp(tab_l, index) {
 	var ligne;
 	var iteration = 0;
 	var i = index;
 	var stop = 0;
 	var t_hsp = new Object();

 	for (i; i < tab_l.length; i++) {
 		ligne = tab_l[i];
 		if (ligne.match(/^>\s*(\S+)/) && stop == 1) {
 			break;
 		} else if (ligne.match(/^>\s*(\S+)/) && stop == 0) {
 			stop = 1;
 		} else if (ligne.match(/^(\s*)?$/)) {
 			continue;
 		}
 		if (ligne.match(/\sScore\s=\s*([0-9]+\.?[0-9]*)/i)) {
 			t_hsp[iteration] = Extract_Hsp(tab_l, i);
 			iteration = iteration + 1;
 		}
 	}
 	return t_hsp;
 }
//*************************************************************************************************************************

 function Extract_Hsp(tab, index) {
 	var hsp = new Object();
 	hsp["Hsp_strand"] = 0;
 	var ligne;
 	var i;
 	var numhsp = 0;
 	var bool = 0;

 	for (i = index; i < tab.length; i++) {
 		ligne = tab[i];
 		if (ligne.match(/\sScore\s=\s*([0-9]+\.?[0-9]*)/i)) {
 			numhsp = numhsp + 1;
 		} else if (ligne.match(/^(\s*)?$/)) {
 			continue;
 		}
 		if (numhsp <= 1) {
 			if (ligne.match(/Score\s*=\s*(\S+)\s*bits\s*(?:\((\d+)\))?,\s*Expect(?:\((\d+\+?)\))?\s*=\s*([^,\s]+)/i)) { //BLAST
 				hsp["Hsp_bit-score"] = RegExp.$1;
 				hsp["Hsp_score"] = RegExp.$2;
 				hsp["Hsp_evalue"] = RegExp.$4;
 				if (hsp["Hsp_evalue"].match(/^(e\S*)/)) {
 					hsp["Hsp_evalue"] = "1" + RegExp.$1;
 				}
 				hsp["Hsp_n"] = RegExp.$3;
 				//??				if(hsp["Hsp_n"]){
 				//??					hsp["Hsp_score"]="";	//deal with BLAST which has no score only bits//??
 				//??				}
 			} else if (ligne.match(/Score\s*=\s*(\S+)\s*\(([\d\.]+)\s*bits\),\s*Expect\s*=\s*([^,\s]+),\s*(?:Sum)?\s*P(?:\(\d+\))?\s*=\s*([^,\s]+)(?:\s*,\s+Group\s*\=\s*(\d+))?/i)) { //wu-BLAST
 				hsp["Hsp_bit-score"] = RegExp.$1;
 				hsp["Hsp_score"] = RegExp.$2;
 				hsp["Hsp_evalue"] = RegExp.$3;
 				hsp["Hsp_pvalue"] = RegExp.$4;
 				hsp["Hsp_group"] = RegExp.$5;
 			} else if (ligne.match(/Identities\s*=\s*(\d+)\s*\/\s*(\d+)\s*[\d\%\(\)]+\s*(?:,\s*Positives\s*=\s*(\d+)\/(\d+)\s*[\d\%\(\)]+\s*)?(?:\,\s*Gaps\s*=\s*(\d+)\/(\d+))?(?:\,\s*Frame\s*=\s*([\+\-][1-3]))?/i)) {
 				hsp["Hsp_identity"] = RegExp.$1;
 				hsp["Hsp_align-len"] = RegExp.$2;
 				if (RegExp.$3) {
 					hsp["Hsp_positive"] = RegExp.$3;
 				} else {
 					hsp["Hsp_positive"] = RegExp.$1;
 				}
 				if (RegExp.$6) {
 					hsp["Hsp_gaps"] = RegExp.$5;
 				}
 				if (RegExp.$7) {
 					hsp["Hsp_frame"] = RegExp.$7;
 					if (RegExp.$7 < 0) {
 						strand = -1
 					}
 					if (RegExp.$7 > 0) {
 						strand = 1
 					}
 				}
 			} else if (ligne.match(/Frame\s*=\s*([\+\-][1-3])\s*(\/\s*([\+\-][1-3]))?/i)) {
 				var Fquery = RegExp.$1;
 				var Fsbjct = RegExp.$2;
 				if (RegExp.$1 && RegExp.$2) {
 					if (Fquery.match(/([\-][1-3])/)) {
 						Squery = -1;
 						strand = -1;
 					}
 					if (Fquery.match(/([\+][1-3])/)) {
 						Squery = 1;
 						strand = 1;
 					}
 					if (Fsbjct.match(/([\-][1-3])/)) {
 						Ssbjct = -1;
 					}
 					if (Fsbjct.match(/([\+][1-3])/)) {
 						Ssbjct = 1;
 					}
 					if (reporttype == "unknow") {
 						reporttype = 'TBLASTX'; // this is for bl2seq only
 					}
 				} else {
 					if (Fquery.match(/([\-][1-3])/)) {
 						Squery = -1;
 					}
 					if (Fquery.match(/([\+][1-3])/)) {
 						Squery = 1;
 					}
 					if (reporttype == "unknow") {
 						reporttype = 'BLASTX';
 					}
 				}
 				if (reporttype) {
 					if (reporttype == 'TBLASTX') {
 						queryframe = RegExp.$1;
 						hitframe = RegExp.$2;
 					}
 					else if (reporttype == "TBLASTN" || reporttype == "PSITBLASTN") {
 						queryframe = 0;
 						hitframe = RegExp.$1;
 					}
 					else if (reporttype == 'BLASTX' || reporttype == 'RPS-BLAST(BLASTP)') {
 						queryframe = RegExp.$1;
 						hitframe = 0;
 						if (reporttype == 'RPS-BLAST(BLASTP)') {
 							reporttype = "RPS-BLAST(BLASTX)";
 						}
 					} 
					else {
 						queryframe = "unknow";
 						hitframe = "unknow";
 					}
 					hsp["Hsp_query-frame"] = queryframe;
 					hsp["Hsp_hit-frame"] = hitframe;
 				}
 			} else if (ligne.match(/Strand\s=\s([a-z]+\s\/\s[a-z]+)/i)) {
 				if (RegExp.$1.match(/(plus\s\/\s[a-z]+)/i)) {
 					hsp["Hsp_strand"] = 1;
 				} else if (RegExp.$1.match(/(minus\s\/\s[a-z]+)/i)) {
 					hsp["Hsp_strand"] = -1;
 				} else {
 					hsp["Hsp_strand"] = RegExp.$1;
 				}
 			} else if (ligne.match(/Links\s*=\s*(\S+)/i)) {
 				hsp["Hsp_links"] = RegExp.$1;
 			} else if (ligne.match(/^((Query|Sbjct):?\s+(\-?\d+)?\s*)(\S+)\s+(\-?\d+)?/) && bool < 2) {
 				var full = RegExp.$1;
 				var type = RegExp.$2;
 				var str = RegExp.$4;
 				var start = RegExp.$3;
 				var end = RegExp.$5;


 				if (type == "Query") {
 					hsp["Hsp_Start_Query"] = start;
 					hsp["Hsp_End_Query"] = end;
 					hsp["Hsp_Str_Query"] = str;
 					var ligne2 = tab[i + 1]; //For the 'homologie' line
 					hsp["Hsp_Str_Homology"] = ligne2;


 					var querylength = str.length;
 					if (hsp["Hsp_Str_Homology"].length != querylength) {
 						var rest = hsp["Hsp_Str_Homology"].length - querylength;
 						hsp["Hsp_Str_Homology"] = hsp["Hsp_Str_Homology"].substring(rest, querylength + rest);
 					}
 				} else if (type == "Sbjct") {
 					hsp["Hsp_Start_Sbjct"] = start;
 					hsp["Hsp_End_Sbjct"] = end;
 					hsp["Hsp_Str_Sbjct"] = str;
 				}
 				bool += 1;
 			} else if (ligne.match(/^((Query|Sbjct):?\s+(\-?\d+)?\s*)(.*)\s+(\-?\d+)?/) && bool == 2) {
 				var full = RegExp.$1;
 				var type = RegExp.$2;
 				var str = RegExp.$4;
				str=str.substring(0,querylength);
 				var end = RegExp.$5;

 				if (type == "Query") {
 					hsp["Hsp_End_Query"] = end;
 					hsp["Hsp_Str_Query"] += str;
 					var ligne2 = tab[i + 1]; 
 					var intermed = ligne2;
 					intermed = intermed.substring(rest, querylength + rest);
 					hsp["Hsp_Str_Homology"] += intermed;
 				} else if (type == "Sbjct") {
 					hsp["Hsp_End_Sbjct"] = end;
 					hsp["Hsp_Str_Sbjct"] += str;
 				}
 			}
 		} else {
 			break;
 		}
 	}
 	if (hsp["Hsp_End_Sbjct"] - hsp["Hsp_Start_Sbjct"] < 0) {
 		var intermediaire = hsp["Hsp_End_Sbjct"];
 		hsp["Hsp_End_Sbjct"] = hsp["Hsp_Start_Sbjct"];
 		hsp["Hsp_Start_Sbjct"] = intermediaire;
 		hsp["the_sframe"] = 0;
 	} else if (hsp["Hsp_End_Sbjct"] - hsp["Hsp_Start_Sbjct"] >= 0) {
 		hsp["the_sframe"] = 1;
 	}
 	if (hsp["Hsp_End_Query"] - hsp["Hsp_Start_Query"] < 0) {
 		var intermediaire = hsp["Hsp_End_Query"];
 		hsp["Hsp_End_Query"] = hsp["Hsp_Start_Query"];
 		hsp["Hsp_Start_Query"] = intermediaire;
 		hsp["the_qframe"] = 0;
 	} else if (hsp["Hsp_End_Query"] - hsp["Hsp_Start_Query"] >= 0) {
 		hsp["the_qframe"] = 1;
 	} else if (hsp["Hsp_frame"] < 0) {
 		var intermediaire = hsp["Hsp_End_Query"];
 		hsp["Hsp_End_Query"] = hsp["Hsp_Start_Query"];
 		hsp["Hsp_Start_Query"] = intermediaire;
 	}
 	if (strand < 0) {
 		hsp["Hsp_strand"] = -1;
 	}
 	if (strand > 0) {
 		hsp["Hsp_strand"] = 1;
 	}
 	if (strand = 0) {
 		hsp["Hsp_strand"] = 0;
 	}
 	return hsp;
 }
//*************************************************************************************************************************

 function Get_Seq_Identifiers(id, ligne) {
 	var acc;
 	if (ligne.match(/(gb|emb|dbj|sp|pdb|bbs|ref|lcl)\|(.*)\|(.*)/)) {
 		var v = RegExp.$2;
 		var elem = v.split(".");
 		acc = elem[0];
 	} else if (ligne.match(/(pir|prf|pat|gnl)\|(.*)\|(\S*)/)) {
 		var v = RegExp.$3;
 		var elem = v.split(".");
 		acc = elem[0];

 	} else {
 		acc = id;
 	}
 	return acc;
 }
//*************************************************************************************************************************

 function Chomp(text) //Need to be create in Javascript
 {
 	if (text != "") {
 		return text.replace(/(\n|\r|\t|\s)*$/, '');
 	}
 }

/*************************************************************************************************************************
**************************************************************************************************************************
**************************************************************************************************************************
**************************************************************************************************************************
*************************************************************************************************************************/

	//////////////////////////////////////////////////////////////////
 	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
 	//:::::::::::::: -Functions to write(HTML)- :::::::::::::::::::://
 	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::://
	//////////////////////////////////////////////////////////////////


//*************************************************************************************************************************


	function Do_Html(){
	 	the_width = window.innerWidth; // Width of  the svg for the graph
	 	var Html_code = "<title>" + blast_hash["BlastOutput_query-def"] + "</title>";
	 	//---------------- CSS ----------------//
	 	Html_code += '<style type="text/css">div.hitdesc{background-color:#EBD5A9;border:1px solid;border-color:grey;};border:1px solid;border-color:grey;}#hit_list{font-size: 0.8em;background-color: #E0E0E0	;}h4{color:blue}all{text-decoration:underline;color:blue}</style>';
	 	//-------------------------------------//

		//---------balise entête---------------//
		var entete=Make_Entete();
		Html_code += '<report id="entete">'+entete+'</report><br>';

	 	//-------------------------------------//

		//---------balise graphique------------//
	 	var nbhit = 10;
	 	var graphique=Make_Graph(nbhit,1); 
		Html_code +='<btn id="Graphic_"><a href="#" onclick=Hide("Graphic")><img src="http://www.aht.li/2075398/moins.png" alt="[-]" height="12" width="12" ></a></btn><br>';
		Html_code += '<report id="Graphic">'+graphique+'</report><br>';
	 	//-------------------------------------//

		//---------balise tabelau de hit-------//
		var Alignement_table=Make_Tab_Hit();
		Html_code +='<btn id="Alignement_table_"><a href="#" onclick=Hide("Alignement_table")><img src="http://www.aht.li/2075398/moins.png" alt="[-]" height="12" width="12"></a></btn>';
		Html_code += '<report id="Alignement_table">'+Alignement_table+'</report>';
	 	//-------------------------------------//

		//---------balise de hits--------------//
	 	conf_href = new Array(); 
		// confirmation table to know if the href of a hit is defined
		var n_hit=0; //to print the first hit

		Html_code +='<br><br><btn id="Alignements_"><a href="#" onclick=Hide("Alignements")><img src="http://www.aht.li/2075398/moins.png" alt="[-]" height="12" width="12"></a></btn>';

		var Alignements=Make_Alignement(n_hit,t_align,1);
		Html_code += '<report id="Alignements">'+Alignements+'</report><br>';
	 	//-------------------------------------//

		//---------balise pied de page---------//
		var pied_de_page=Make_Pied_De_Page();
		Html_code += '<report id="pied_de_page">'+pied_de_page+'</report>';
	 	//-------------------------------------//

	 	document.body.innerHTML = Html_code; 
	}

//*************************************************************************************************************************

function Make_Entete(){

	var Html_code="";

	//---------------- Links ----------------//
 	Html_code += '<b id="top">Query</b>= ' + '<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=search&amp;term=' + blast_hash["BlastOutput_query-def"] + '" target="_blank" title="NCBI link">' + blast_hash["BlastOutput_query-def"] + "</a> (" + blast_hash["BlastOutput_query-len"] + " letters) [" + blast_hash["BlastOutput_program"] + "]";
 	Html_code += "<br><b>Database:</b> " + blast_hash["Parameters_full_dbpath"] + "; " + blast_hash["BlastOutput_db-len"] + " sequences; " + blast_hash["BlastOutput_db-let"] + " total letters<br>";
 	Html_code += '<ul><li><a href="#hit_list">Alignments list</a></li><li><a href="#alignments">Alignments</a></li><li><a href="#parameters">Search parameters & Bottom</a></li><li><a href="http://www.aht.li/2079424/The_help.html" target="_blank">Help<img src="http://www.aht.li/2079435/inte.svg" alt="?"  height="12" width="15"></a></li></ul>';




 	//-------------------------------------//
	//-------------- Options ---------------//

	Html_code +='<btn id="Visual_Parameters_"><a href="#" onclick=Hide("Visual_Parameters")><img src="http://www.aht.li/2075398/moins.png" alt="[-]" height="12" width="12"></a></btn>';
 	Html_code += '<div style="background-color:#D1C2B2;border:1px solid;width:500px;"  id="Visual_Parameters" class="Visual_Parameters" ><BLOCKQUOTE> <h3>Visual Parameters</h3></BLOCKQUOTE><ul>';

	Html_code += '<li>Change the alignement size: <input type="text" name="input1" style="width: 30px"  onchange="Make_Alignement(0,this.value)" id="taille_alignements"/>[ <a href="#" onclick="Make_Alignement(0,60)">60</a> default]</li>';
 	Html_code += '<li>Display <input type="text" name="input2" style="width: 30px"  onchange="Do_Some_Hits(this.value,' + t_align + ')" id="nombre_hits"/>HITs (<a href="#" onclick="Do_All_Hits()">' + (Object.keys(blast_hash["BlastOutput_tabHit"])).length + ' max</a>)</li>';
 	Html_code += '<li>Display <input type="text" name="input3" style="width: 30px"  onchange="Make_Graph(this.value)" id="nombre_hits_graph"/>HITs on the graph (<a href="#" onclick="Make_Graph(' + (Object.keys(blast_hash["BlastOutput_tabHit"])).length + ')">' + (Object.keys(blast_hash["BlastOutput_tabHit"])).length + ' max</a>)</li>';

	Html_code += '<li>Display <input type="text" name="input4" style="width: 30px"  onchange="nbhsp=this.value;Make_Graph(10)" id="nombre_hits_graph"/>HSPs by HITs on the graph (<a href="#" onclick="nbhsp=10000;Make_Graph(10)">max</a>)</li>';

	Html_code += '<li>Minimum score for the display of HSPs <input type="text" name="input5" style="width: 30px"  onchange="scorehsp=parseInt(this.value);Make_Graph(10)" id="score_hsp_graph"/>[ <a href="#" onclick="scorehsp=0;Make_Graph(10)">0</a> default]</li>';


	Html_code += '</ul></div><br>';
 	//--------------------------------------//


	return Html_code;

}


//*************************************************************************************************************************
function Make_Graph(nbhit,where){
	var Html_code="";

 	var the_width2 = the_width - 100;
 	var fact = the_width2 / blast_hash["BlastOutput_query-len"];
 	var pos = 40;

 	//-------------- Scale -------------------//
 	Html_code += Do_Scale(fact);
 	//----------------------------------------//

 	//---------------- HITs -----------------//
 	for (var j = 0; j < (Object.keys(blast_hash["BlastOutput_tabHit"])).length && j < nbhit; j++) { 
 		var pos = pos + 30;

 		Html_code += '<t  id="graph_hit' + j + '">'
 		var hit = Draw_Hit(j, fact, 1);
 		Html_code += hit;
 		Html_code += '</t>'
 	}
 	//--------------------------------------//

 	//------------ Score Scale -------------//
Html_code += Do_Color_Scale();
 	//--------------------------------------//



 	if (where == 1) { //permet de savoir si l'utilisation du graph est la première ou non
		return Html_code;
 	} else { //sinon l'appel est effectué par le paramètre et on change ainsi la balise graph
 		document.getElementById('Graphic').innerHTML = Html_code;
 	}
}

//*************************************************************************************************************************
function Make_Tab_Hit(){
	var Html_code="";

	//-------------- Tab of Hit ---------------//
 	Html_code += '<br><table border=0 id="hit_list" CELLPADDING=3><tr><th>N&deg;</th><th>Sequences producing significant alignments:</th><th>Score<br>(bits)</th><th>E<br>value</th><th>nhsp</th></tr>';
 	for (var j = 0; j < (Object.keys(blast_hash["BlastOutput_tabHit"])).length; j++) {
 		var n_hit = j + 1;

 		Html_code += '<tr><td><a name="list_' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '"></a><a href="#hit_'+j+'" onclick="Go_To_Hit(' + j + ')">' + " -" + n_hit + "- </a></td><td>" + '<a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=' + 'nucleotide' + '&amp;cmd=search&amp;term=' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '" target="_blank" title="NCBI link">' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '</a> ' + blast_hash["BlastOutput_tabHit"][j]["Hit_def"] + "</td> <td ALIGN=RIGHT>" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][0]["Hsp_score"] + '</td><td ALIGN=RIGHT>' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][0]["Hsp_evalue"] + '</td><td ALIGN=RIGHT>' + Object.keys(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"]).length + '</td></tr>';
 	}
 	Html_code += '</table>';
	//-----------------------------------------//

	return Html_code;
}

//*************************************************************************************************************************
function Make_Pied_De_Page(){
	var Html_code="";
 
	//-------------- Search Parameters ---------------//
if(blast_hash["BlastOutput_algorithm-reference"]){
	if (blast_hash["BlastOutput_algorithm-reference"].match(/"(.*)"/i)) {
 			var ref= RegExp.$1;
			ref=ref.replace(/\s/g, '+');
	}
}
 	Html_code += '<hr><h2 id="parameters">Search Parameters</h2>' + "<table border=1><tr><th>Parameter</th><th>Value</th></tr>" + '<tr><td>allowgaps</td><td>' + blast_hash["Parameters_allowgaps"] + "</td></tr>" + "<tr><td>expect</td><td>" + blast_hash["Parameters_expect"] + "</td></tr>" + "<tr><td>gapext</td><td>" + blast_hash["Parameters_gap-extend"] + "</td></tr>" + "<tr><td>gapopen</td><td>" + blast_hash["Parameters_gap-open"] + "</td></tr>" + "<tr><td>matrix</td><td>" + blast_hash["Parameters_matrix"] + "</td></tr>" + "</table><br><hr>";
 	Html_code += '<a href="#top"><img src="http://www.aht.li/2076124/fleche.svg" alt="/\\" height="28" width="28" ></a><br>';
 	Html_code += '<br><a href=http://www.ncbi.nlm.nih.gov/pubmed/?term='+ref+' target="_blank" title="PubMed link">' + blast_hash["BlastOutput_algorithm-reference"] + '</a><br>';


Html_code += '<br><a href="#" onclick="Recup_html();">Capture</a>';


var aide="On click it show the HTML code of the current page.<br>Copy and paste the code in a new document and save it on the HTML format. <br>The document will be as you capture. <br>If you open your HTML document, only links to the ncbi will work (unless you open all the hits) and the graphic will not be interactive.<br>/!\\ Care for the big reports./!\\";

Html_code += '<a href="#parameters" onclick=document.getElementById("help_capture").style.display="block";><img src="http://www.aht.li/2079435/inte.svg" alt="?"  height="12" width="15"></a><br><div id="help_capture" style="background-color:#FFDEAE">'+aide+'</div>';

 	//-----------------------------------------------//
	//-------------- Information ---------------//
 	Html_code +='<h5>Blast parsing and representation by <a href="http://www.aht.li/2062427/JETBlastReport.png" target="_blank">JETBlastReport</a> [Javascript Enhancer & Transformer of Blast Reports]<br>Version of  <font color="purple" >benoit.piegu&nbsp;@&nbsp;tours.inra.fr</FONT> and <font color="purple" > valentin.marcon&nbsp;@&nbsp;etu.udamail.fr</font> (2013/06/27)</h5>';
 	//------------------------------------------//

	return Html_code;
}

//*************************************************************************************************************************
function Make_Alignement(n_hit,t_align2,when) {

	t_align=t_align2;
 	var Html_code = '<a name="hit_'+n_hit+'"></a>';
 	var j = n_hit;
 	var ind = n_hit+1;

 	conf_href[ind] = "ok"; //confirme que le lien vers le hit est definit

 	Html_code += '<a name="alignments"></a>' + '<br><div id="' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '" class="hitdesc" >' + '<a href="#top">^</a><br>' + '<a href="#list_' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '">' + ind + '</a>' + '><b><a href="http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&amp;cmd=search&amp;term=' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '" target="_blank" title="NCBI link">' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '</a></b> ' + blast_hash["BlastOutput_tabHit"][j]["Hit_def"] + "<br>" + '<strong>Length</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["Hit_len"] + "</div><br>";
 	var l = 0;

 	//// HSP
 	for (var i = 0; i < Object.keys(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"]).length; i++) {
 		var ql_aln = Math.round(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] / blast_hash["BlastOutput_query-len"] * 100);
 		var sl_aln = Math.round(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] / blast_hash["BlastOutput_db-let"] * 100);
//***
 		Html_code += '<a name="hit_'+n_hit+'_'+l+'"></a><br><div id="hsp_'+j+'_'+i+'_" ><strong><a href="#ghsp_'+n_hit+'">G</a>/<a href="#list_' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + '">T</a>-Length</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] + "(%ql_aln: " + ql_aln + "%, %sl_aln: " + sl_aln + "%), <strong>Score</strong> = " + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_bit-score"] + " bits(" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_score"] + "), <strong>Expect = </strong>" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_evalue"] + "<br>" + '<strong>Identities</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_identity"] + "/" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] + " (" + Math.round(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_identity"] / blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] * 100) + "%), ";
 		if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_hit-frame"] != "unknow" && typeof (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_hit-frame"]) != 'undefined') {
 			Html_code += '<strong>Frame</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_hit-frame"] + ", ";
 		}
 		if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_gaps"] != "unknow" && typeof (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_gaps"]) != 'undefined') {
 			Html_code += '<strong>Gaps</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_gaps"] + "/" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] + " (" + Math.round(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_gaps"] / blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] * 100) + "%), ";
 		}
 		if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_positive"] != "unknow" && typeof (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_positive"]) != 'undefined' && blast_hash["BlastOutput_program"] != "BLASTN") {
 			Html_code += '<strong>Positives</strong> = ' + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_positive"] + "/" + blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] + " (" + Math.round(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_positive"] / blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_align-len"] * 100) + "%)";
 		}


 		//// SEQ
 		var x = 0;
 		var y = 0;
 		var z = 0;

 		if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Sbjct"].length >= blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Query"].length) {
 			var espace = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Sbjct"].length;
 		} else {
 			var espace = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Query"].length;
 		}
 		Html_code += "<pre>";
 		var coefQ= 1;
 		var coefS = 1;

 		if (blast_hash["BlastOutput_program"] == "TBLASTX") {
 			coefQ = 3;
 			coefS = 3;
 		} else if (blast_hash["BlastOutput_program"] == "BLASTX") {
 			coefQ = 3;
 		} else if (blast_hash["BlastOutput_program"] == "TBLASTN") {
 			coefS = 3;
 		}
		var length=blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Query"]-blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Start_Query"];
 		while (y * coefQ < length) {

 			var str_query = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Str_Query"].substr(y, t_align);
 			var str_homol = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Str_Homology"].substr(y, t_align);
 			var str_sbjct = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Str_Sbjct"].substr(y, t_align);
 			var l_query = str_query.replace(/-/g, '').length * coefQ - 1;
 			var l_sbjct = str_sbjct.replace(/-/g, '').length * coefS - 1;
 			if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["the_qframe"] > 0) {
 				var deb_query = parseInt(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Start_Query"]) + parseInt(x);
 				var fin_query = deb_query + l_query;
 			} else if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["the_qframe"] == 0) {
 				var deb_query = parseInt(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Query"]) - parseInt(x);
 				var fin_query = deb_query - l_query;
 			}
 			if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["the_sframe"] > 0) {
 				var deb_sbjct = parseInt(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_Start_Sbjct"]) + parseInt(z);
 				var fin_sbjct = deb_sbjct + l_sbjct;
 			} else if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["the_sframe"] == 0) {
 				var deb_sbjct = parseInt(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][l]["Hsp_End_Sbjct"]) - parseInt(z);
 				var fin_sbjct = deb_sbjct - l_sbjct;
 			}
 			Html_code += "<br>" + 'Query: ' + deb_query + " ";
 			var espace2 = deb_query.toString().length;
 			var esp_query = Write_Caracter(" ", espace - espace2);
 			Html_code += esp_query + str_query + " " + fin_query + "<br>";
 			var esp_homol = Write_Caracter(" ", espace + 8);
 			Html_code += esp_homol + str_homol + "<br>" + "Sbjct: ";
 			Html_code += deb_sbjct + " ";
 			var espace3 = deb_sbjct.toString().length;
 			var esp_sbjct = Write_Caracter(" ", espace - espace3);
 			Html_code += esp_sbjct + str_sbjct + " " + fin_sbjct + "<br>";
 			x = x + l_query + 1;
 			y = y + parseInt(t_align);
 			z = z + l_sbjct + 1;
 		}
 		l = l + 1;
		Html_code += '</div>';
 		Html_code += "</pre>";
 	}

 	if (j + 1 < (Object.keys(blast_hash["BlastOutput_tabHit"])).length && when!=3) {
 		Html_code += '<br><hit id="myHit' + j + '"><h4 onclick="Do_Next_Hit(' + j + ')" >>>[' + blast_hash["BlastOutput_tabHit"][j+1]["Hit_accession"] + ']</h4>';
 		Html_code += '<h4 href="#" onclick="Do_All_Hits()">>>[All Hits]</h4></hit>';
 	}// To have the next hit
 	Html_code += '<hit id="Bonus_hit"></hit>'; // To print an extra hit, without print all of them
 	




	if(when==1 || when==3){
 		return Html_code;
	}
	else{
 		document.getElementById('Alignements').innerHTML =Html_code;
	}

}

//*************************************************************************************************************************

 function Do_Next_Hit( n_hit) {
 	document.getElementById('Bonus_hit').innerHTML = "";
 	var j = n_hit + 1;
 	var Html_code = Make_Alignement( j, t_align ,1);

 	document.getElementById('myHit' + n_hit).innerHTML = Html_code;

 }
//*************************************************************************************************************************

 function Do_Some_Hits(nb_hit) {
	
	if(nb_hit>1){
 		for (var j = 0; j < nb_hit - 1; j++) {
 			if (nb_hit <= (Object.keys(blast_hash["BlastOutput_tabHit"])).length) {
 				Do_Next_Hit( j);
 			}
 		}
	} 
	else{
		Do_Html();
	}
 }
//*************************************************************************************************************************

 function Go_To_Hit(j) {
	n=j+1;
 	if (conf_href[n] != "ok") {
 		// Si le hit n'est pas présent il est inséré dans une balise bonus
 		var Html_code = '<hr>';
 		Html_code += Make_Alignement( j,t_align,3);
 		Html_code += '<br>';
 		conf_href[n] = "0"; //cause its just a bonus
 		document.getElementById('Bonus_hit').innerHTML = Html_code;
 	}
 }
//*************************************************************************************************************************

function Go_To_Hsp(j,i){

		Go_To_Hit(j);
		var id='hsp_'+j+'_'+i+'_';
		document.getElementById(id).style.backgroundColor="#94E994";
		var interv=setInterval(function(){document.getElementById(id).style.backgroundColor="#FFFFFF";window.clearInterval(interv)},1500);

}
//*************************************************************************************************************************

 function Do_All_Hits() {
 	for (var j = 0; j < (Object.keys(blast_hash["BlastOutput_tabHit"])).length; j++) {
 		Do_Next_Hit( j);
 	}
 }
//*************************************************************************************************************************

 function Write_Caracter(caracter, n) {
 	var string = "";
 	for (var i = 0; i < n; i++) {
 		string += caracter;
 	}
 	return string;
 }

/************************************************************************************************************************
*************************************************************************************************************************
*************************************************************************************************************************/

	//::::::::::::::::::::::::::::::::::::::::::::::::://
	//::::::::::Function to make the graph::::::::::::://
	//::::::::::::::::::::::::::::::::::::::::::::::::://



 function Draw_Hit(j, fact, where) {
 	var the_width2 = the_width - 100;
 	var code = '<br><a name="ghsp_'+j+'"></a><svg width="' + the_width + '" height="35" id="svg_graph_hit' + j + '">';
 	var posy = 12;
 	var mult_hsp = "";
	for(var i=0;i<tab_cord.length; i++){
 					code += '<line x1="' +tab_cord[i] + '" y1="0" x2="' + tab_cord[i]  + '" y2="35" style="fill:none;stroke:#CC9900;stroke-width:0.5" />';
	}
 	if (Object.keys(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"]).length > 1) {
		var posx=3+(9*blast_hash["BlastOutput_tabHit"][j]["Hit_accession"].length);

		code += '<a xlink:href="#" onclick="Draw_Hsp_Of_Hit(' + j + ',' + fact + ')" >	<image xlink:href="http://www.aht.li/2076119/fleche2.svg" x="4" y="2" height="10px" width="10px"/></a>'; 



 	} 
 	code += '<a xlink:href="#hit_'+j+'" onclick="Go_To_Hit(' + j + ')" ><text x="20" y="' + posy + '" fill="green"   >' + blast_hash["BlastOutput_tabHit"][j]["Hit_accession"] + mult_hsp + '</text></a>' +"\n"; 

 	var posy = 19;
 	code += ' <line x1="0" y1="' + posy + '" x2="' + the_width2 + '" y2="' + posy + '" style="fill:none;stroke:black;stroke-dasharray:3,2;stroke-width:0.2" />' +"\n";


	if (nbhsp>Object.keys(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"]).length){
		var nbhsp2=Object.keys(blast_hash["BlastOutput_tabHit"][j]["tab_hsp"]).length;
	}
	else{
		var nbhsp2=nbhsp;
	}
 	for (var i = 0; i < nbhsp2; i++) { //HSP
 		var score = blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][i]["Hsp_score"];
		if(score<scorehsp){
 			continue;
		}
 		var col = Get_Score_Color(score);
 		var width = (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][i]["Hsp_End_Query"] - blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][i]["Hsp_Start_Query"]) * fact
 		var deb = (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][i]["Hsp_Start_Query"] * fact);
 		var y = 15;
 		var y2 = y+ 4;
 		var y3 = y + 8;
 		var x = deb;
 		var x3 = x + width;

 		if (width >= 5) {
 			var x4 = x + 5;
 			var x2 = x + width - 5;
 		} else {
 			var x4 = x3;
			var x2=x;
 		}

	code +='<a xlink:href="#hit_'+j+'_'+i+'"  onclick="Go_To_Hsp(' + j + ','+i+')"  >';
 		if (blast_hash["BlastOutput_tabHit"][j]["tab_hsp"][i]["the_sframe"] == 1) {
 			code += '<polyline points="' + x + ',' + y + ' ' + x2 + ',' + y + ' ' + x3 + ',' + y2 + ' ' + x2 + ',' + y3 + ' ' + x + ',' + y3 + ' ' + x + ',' + y + '"  fill="' + col + '" style="stroke:black;stroke-width:0.5"  />';
 		} else {
 			code += '<polyline  points="' + x + ',' + y2 + ' ' + x4 + ',' + y + ' ' + x3 + ',' + y + ' ' + x3 + ',' + y3 + ' ' + x4 + ',' + y3 + ' ' + x + ',' + y2 + '"  fill="' + col + '" style="stroke:black;stroke-width:0.5" />' +"\n";
 		}
		code+='</a>';

 	}
 	if (where == 1) {
 		code += '</svg>';
 		return code;
 	} else {
 		var id = 'graph_hit' + j;
 		document.getElementById(id).innerHTML = code;
 		document.getElementById(id).target = "_blank";
 	}
 }
//*************************************************************************************************************************

 function Draw_Hsp_Of_Hit(n, fact) {

 	var taille = 26;

	if (nbhsp>Object.keys(blast_hash["BlastOutput_tabHit"][n]["tab_hsp"]).length){
		var nbhsp2=Object.keys(blast_hash["BlastOutput_tabHit"][n]["tab_hsp"]).length;
	}
	else{
		var nbhsp2=nbhsp;
	}

 	taille += 9 * nbhsp2;
 	var posy = 12;
 	var the_width2 = the_width - 100
 	var mult_hsp = "";
 	var code = '<a name="ghsp_'+n+'"></a><svg width="' + the_width + '" height="' + taille + '" >';

	for(var i=0;i<tab_cord.length; i++){
 					code += '<line x1="' +tab_cord[i] + '" y1="0" x2="' + tab_cord[i]  + '" y2="'+taille+'" style="fill:none;stroke:#CC9900;stroke-width:0.5" />';
	}

 	if (Object.keys(blast_hash["BlastOutput_tabHit"][n]["tab_hsp"]).length > 1) {
		var posx=3+(9*blast_hash["BlastOutput_tabHit"][n]["Hit_accession"].length);

		code += '<a xlink:href="#" onclick="Draw_Hit(' + n + ',' + fact + ',0)" >	<image xlink:href="http://www.aht.li/2076120/fleche3.svg" x="4" y="2" height="10px" width="10px"/></a>'; 
 	} 
 	code += '<a xlink:href="#hit_'+n+'" onclick="Go_To_Hsp(' + n + ','+i+')" ><text x="20" y="' + posy + '" fill="green"   >' + blast_hash["BlastOutput_tabHit"][n]["Hit_accession"] + mult_hsp + '</text></a>'; 

	posy=6;

 	for (var i = 0; i < nbhsp2; i++) { //HSP
 		var score = blast_hash["BlastOutput_tabHit"][n]["tab_hsp"][i]["Hsp_score"];
		if(score<scorehsp){
 			continue;
		}
 		var col = Get_Score_Color(score);
 		var width = (blast_hash["BlastOutput_tabHit"][n]["tab_hsp"][i]["Hsp_End_Query"] - blast_hash["BlastOutput_tabHit"][n]["tab_hsp"][i]["Hsp_Start_Query"]) * fact;
 		var deb = (blast_hash["BlastOutput_tabHit"][n]["tab_hsp"][i]["Hsp_Start_Query"] * fact);
 		posy = posy + 13;
 		code += ' <line x1="0" y1="' + posy + '" x2="' + the_width2 + '" y2="' + posy + '" style="fill:none;stroke:black;stroke-dasharray:3,2;stroke-width:0.2" />';
		posy = posy -4;
 		var y = posy;
 		var y2 = posy +4 ;
 		var y3 = posy +8;
 		var x = deb;
 		var x3 = x + width;

 		if (width >= 5) {
 			var x4 = x + 5;
 			var x2 = x + width - 5;
 		} else {
 			var x4 = x3;
			var x2=x;
 		}
	code +='<a xlink:href="#hit_'+n+'_'+i+'" onclick="Go_To_Hsp(' + n + ','+i+')" >';
 		if (blast_hash["BlastOutput_tabHit"][n]["tab_hsp"][i]["the_sframe"] == 1) {
 			code += '<polyline  points="' + x + ',' + y + ' ' + x2 + ',' + y + ' ' + x3 + ',' + y2 + ' ' + x2 + ',' + y3 + ' ' + x + ',' + y3 + ' ' + x + ',' + y + '"  fill="' + col + '" style="stroke:black;stroke-width:0.5"  />';
 		} else {
 			code += '<polyline  points="' + x + ',' + y2 + ' ' + x4 + ',' + y + ' ' + x3 + ',' + y + ' ' + x3 + ',' + y3 + ' ' + x4 + ',' + y3 + ' ' + x + ',' + y2 + '"  fill="' + col + '" style="stroke:black;stroke-width:0.5" />';
 		}
		code+='</a>';
 	}
 	code += '</svg>';
 	var id = 'graph_hit' + n;
 	document.getElementById('graph_hit' + n).innerHTML = code;
 	document.getElementById('graph_hit' + n).target = "_blank";
 }

//*************************************************************************************************************************

 function Get_Score_Color(score) {
 	var color = "";

 	if (score >= 0) {
 		if (score <= 20) {
 			color = "#FFFFFF";
 		} else if (score <= 40) {
 			color = "#FFDFDF";
 		} else if (score <= 60) {
 			color = "#FFCACA";
 		} else if (score <= 80) {
 			color = "#FFBDBD";
 		} else if (score <= 100) {
 			color = "#FFACAC";
 		} else if (score <= 120) {
 			color = "#FF9797";
 		} else if (score <= 140) {
 			color = "#FF7D7D";
 		} else if (score <= 160) {
 			color = "#FF5C5C";
 		} else if (score <= 180) {
 			color = "#FF3333";
 		} else if (score > 180) {
 			color = "#FF0000";
 		}
 	}
 	return color;
 }

//*************************************************************************************************************************

 function Do_Color_Scale() {
 	var the_width2 = the_width - 100;
 	var taille = 15;
 	var taille2 = taille + 24;
 	var s = 0;
 	var code = '<br><svg width="' + the_width + '" height="45" id="colorScale">';
	for(var i=0;i<tab_cord.length; i++){
 					code += '<line x1="' +tab_cord[i] + '" y1="0" x2="' + tab_cord[i]  + '" y2="45" style="fill:none;stroke:#CC9900;stroke-width:0.5" />';
	}
 	var tab_col = ["#FFFFFF", "#FFDFDF", "#FFCACA", "#FFBDBD", "#FFACAC", "#FF9797", "#FF7D7D", "#FF5C5C", "#FF3333", "#FF0000"];
 	for (var i = 0; i < 10; i++) {
 		code += '<rect width=' + the_width2 / 10 + ' height="10" x=' + i * the_width2 / 10 + ' y="' + taille + '" fill="' + tab_col[i] + '" style="stroke:black;stroke-width:1"/>';
 		code += '<text x="' + i * the_width2 / 10 + '" y="' + taille2 + '" fill="#000000">s>' + s + ' </text>';
 		s = s + 20;
 	}
 	code += '<svg><br>';
 	return code;
 }
//*************************************************************************************************************************

 function Do_Scale(fact) {
 	var tab_bornes = [1, 10, 100, 200, 1000, 2000, 10000, 20000, 100000, 200000, 400000, 800000, 160000];
 	var code = "";
	var z=0;
	tab_cord=new Array();

	if(blast_hash["BlastOutput_query-len"]>160000){ //if the query len is bigger than 160 000
		for(i=160000;i<blast_hash["BlastOutput_query-len"];i=i*2){
			j=i*2;
			var tab_bornes = [i,j];
		}
	}

 	for (var i = 0; i < (tab_bornes.length) - 1; i++) {
 		var a = tab_bornes[i];
 		var b = tab_bornes[i + 1];
 		if (blast_hash["BlastOutput_query-len"] <= b && blast_hash["BlastOutput_query-len"] > a) {
 			var test = 0;
 			var the_width2 = the_width - 100;
 			var code = '<svg width="' + the_width + '" height="40" id="colorScale">';
 			code += ' <line x1="1" y1="10" x2="' + the_width2 + '" y2="10" style="fill:none;stroke:blue;stroke-width:0.8" />';
 			var c = (b + a) / 2;
 			var a2 = a / 10;
 			var b2 = b / 10;

 			for (var i = 0; i < blast_hash["BlastOutput_query-len"]; i = i + a2) {
 				var espace = 1 + i * fact;
 				var espace2 = espace - 12;

 				if (test == b2) {
 					test = 0;
 					code += '<text x="' + espace2 + '" y="35" fill="blue">' + i + '</text>';
 					code += '<line x1="' + espace + '" y1="0" x2="' + espace + '" y2="20" style="fill:none;stroke:blue;stroke-width:0.8" />';
 				}
 				test = test + a2;
 				code += '<line x1="' + espace + '" y1="3" x2="' + espace + '" y2="17" style="fill:none;stroke:blue;stroke-width:0.8" />';
 				code += '<line x1="' + espace + '" y1="0" x2="' + espace + '" y2="40" style="fill:none;stroke:#CC9900;stroke-width:0.5" />';
				tab_cord[z++]=espace;
 			}

 			if (blast_hash["BlastOutput_query-len"] <= c) {
 				var b3 = b / 100;
 				for (var i = 0; i < blast_hash["BlastOutput_query-len"]; i = i + b3) {
 					var espace = 1 + i * fact;
 					var espace2 = espace - 12;
 					code += '<line x1="' + espace + '" y1="0" x2="' + espace + '" y2="40" style="fill:none;stroke:#CC9900;stroke-width:0.5" />';
					tab_cord[z++]=espace;
 					code += '<line x1="' + espace + '" y1="6" x2="' + espace + '" y2="14" style="fill:none;stroke:blue;stroke-width:0.8" />';
 				}
 			}
 			code += '</svg>';
 			return code;
 		}
 	}
	

 }


//*************************************************************************************************************************

function Display(obj){
        document.getElementById(obj).style.display = "block";
				var obj2=obj+"_";
				var btnm='<a href="#" onclick=Hide("'+obj+'")><img src="http://www.aht.li/2075398/moins.png" alt="[-]" height="12" width="12"></a>';
				document.getElementById(obj2).innerHTML=btnm
}
 
function Hide(obj){
        document.getElementById(obj).style.display = "none";
				var obj2=obj+"_";
				var btnp='<a href="#"  onclick=Display("'+obj+'") ><img src="http://www.aht.li/2075430/plus.png" alt="[-]" height="12" width="12"> <b>'+obj+'</b></a>';
				document.getElementById(obj2).innerHTML=btnp
}

//*************************************************************************************************************************

function Recup_html(){
var all=" <xmp>"; //html is not interpreted
all+=document.documentElement.innerHTML;
all+="</xmp>";
document.documentElement.innerHTML=all;
}

//*************************************************************************************************************************

