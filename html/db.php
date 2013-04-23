<?php
	if (isset($_POST["download"]))
	{
		header("Content-Type: text/plain");
		header("Content-Disposition: attachment; filename=\"" . $_POST["qty"] . ".txt\"");
		echo $_POST["hexoutput"];
		exit();
	}
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
 
<head>
	<title>Hex scattering database</title>
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
	
	<link rel="icon" type="image/gif" href="hexe-small.gif" />
	
	<style type="text/css">
		<!-- @import "style.css"; -->
	</style>

	<script>
		function jsAvailClick()
		{
			if (document.getElementById("avail-more-1").style.display == "none")
			{
				// change the title text
				document.getElementById("avail-link").innerHTML = "&#x25B2; Available data &#x25B2;"
				
				// display info
				document.getElementById("avail-more-1").style.display = "block";
				document.getElementById("avail-more-2").style.display = "block";
				document.getElementById("avail-more-3").style.display = "block";
			}
			else
			{
				// change the title text
				document.getElementById("avail-link").innerHTML = "&#x25BC; Available data &#x25BC;"
				
				// hide info
				document.getElementById("avail-more-1").style.display = "none";
				document.getElementById("avail-more-2").style.display = "none";
				document.getElementById("avail-more-3").style.display = "none";
			}
		}
		function jsContextClick()
		{
			if (document.getElementById("context-more-1").style.display == "none")
			{
				// change the title text
				document.getElementById("context-link").innerHTML = "&#x25B2; Data context &#x25B2;"

				// display info
				document.getElementById("context-more-1").style.display = "block";
			}
			else
			{
				// change the title text
				document.getElementById("context-link").innerHTML = "&#x25BC; Data context &#x25BC;"

				// display info
				document.getElementById("context-more-1").style.display = "none";
			}
		}
	</script>

	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
			extensions: ["tex2jax.js"],
			jax: ["input/TeX","output/HTML-CSS"],
// 			tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
		});
	</script>

	<script src="http://www.mathjax.org/mathjax/MathJax.js"></script>
</head>
 
<body>

	<div class = "grid">

	<!-- header -->

	<center>
	<table border = "0">
		<tr><td>
			<a href = "index.html"><img src = "hexe.gif" border = "0" alt = "logo"/></a>
		</td><td>
			<div class = "nadpis">Hex</div>
			
			<div class = "podnadpis">scattering database</div>
		</td></tr>
	</table>
	</center>

	<table width = "100%" style = "table-layout:fixed;"><tr><td valign = "top" width = "50%">

	<div class = "sekce">Input:</div>

	<form name = "data" action = "db.php" method = "post" style = "margin: 10px">
		<div class = "text">Choose what to compute:</div>
		<center>
			<select name = "qty" onchange="this.form.submit()">
				<option value = "ccs" <?php
					if (!isset($_POST["qty"]) or $_POST["qty"] == "ccs")
						echo "selected = \"selected\"";
					?> >complete cross section</option>
				<option value = "colls" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "colls")
						echo "selected = \"selected\"";
					?> >collision strength</option>
				<option value = "dcs" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "dcs")
						echo "selected = \"selected\"";
					?> >differential cross section</option>
				<option value = "xcs" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "xcs")
						echo "selected = \"selected\"";
					?> >extrapolated cross section</option>
				<option value = "ics" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "ics")
						echo "selected = \"selected\"";
					?> >integral cross section</option>
				<option value = "momtf" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "momtf")
						echo "selected = \"selected\"";
					?> >momentum transfer</option>
				<option value = "scatamp" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "scatamp")
						echo "selected = \"selected\"";
					?> >scatering amplitude</option>
				<option value = "tcs" <?php
					if (isset($_POST["qty"]) and $_POST["qty"] == "tcs")
						echo "selected = \"selected\"";
					?> >total cross section</option>
			</select>
		</center>

		<div class = "text">Choose units:</div>
		<center>
			<span class = "text">Energy: <select name = "Eunits">
				<option value = "Ry" <?php
					if (!isset($_POST["Eunits"]) or $_POST["Eunits"] == "Ry")
						echo "selected = \"selected\"";
					?> >Rydberg</option>
				<option value = "a.u." <?php
					if (isset($_POST["Eunits"]) and $_POST["Eunits"] == "a.u.")
						echo "selected = \"selected\"";
					?> >Hartree (a.u.)</option>
				<option value = "eV" <?php
					if (isset($_POST["Eunits"]) and $_POST["Eunits"] == "eV")
						echo "selected = \"selected\"";
					?> >eV</option>
			</select> ,  Output: <select name = "Tunits">
				<option value = "a.u." <?php
					if (!isset($_POST["Tunits"]) or $_POST["Tunits"] == "a.u.")
						echo "selected = \"selected\"";
					?> >a.u.</option>
				<option value = "cgs" <?php
					if (isset($_POST["Tunits"]) and $_POST["Tunits"] == "cgs")
						echo "selected = \"selected\"";
					?> >cgs</option>
			</select></span>
		</center>

		<div class = "text">Set initial atomic state(s):</div>
		<center>
			\(n_i\) = <input type = "text" name = "ni" size = "3" value = "<?php echo (isset($_POST["ni"]) ? $_POST["ni"] : 1); ?>"/>
			\(l_i\) = <input type = "text" name = "li" size = "3" value = "<?php echo (isset($_POST["li"]) ? $_POST["li"] : 0); ?>"/>
			\(m_i\) = <input type = "text" name = "mi" size = "3" value = "<?php echo (isset($_POST["mi"]) ? $_POST["mi"] : 0); ?>"/>
		</center>
		
		<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] != "tcs")
			{
				printf("\t\t<div class = \"text\">Set final atomic state(s):</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(n_f\\) = <input type = \"text\" name = \"nf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["nf"]) ? $_POST["nf"] : "1");
				printf("\t\t\t\\(l_f\\) = <input type = \"text\" name = \"lf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["lf"]) ? $_POST["lf"] : "0");
				printf("\t\t\t\\(m_f\\) = <input type = \"text\" name = \"mf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["mf"]) ? $_POST["mf"] : "0");
				printf("\t\t</center>\n");
			}
		?>
		
		<?php
			if (isset($_POST["qty"]) and $_POST["qty"] == "scatamp") //~ scatamp
			{
				printf("\t\t<div class = \"text\">Set global quantum numbers:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(E\\) = <input type = \"text\" name = \"E\" size = \"3\" value = \"%s\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
				printf("\t\t\t\\(S\\) = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
				printf("\t\t</center>\n");
			}
			else if (!isset($_POST["qty"]) or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "tcs")
			{
				// do nothing
			}
			else
			{
				printf("\t\t<div class = \"text\" title = \"'L' is the total angular momentum, 'S' is the total spin.\">Set global quantum numbers:</div>\n");
				printf("\t\t<center>\n");
				switch ($_POST["qty"])
				{
					case "scatamp":
					case "dcs":
						printf("\t\t\t\\(E\\) = <input type = \"text\" name = \"E\" size = \"3\" value = \"%s\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
						printf("\t\t\t\\(S\\) = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "momtf":
						printf("\t\t\t\\(S\\) = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "ics":
					case "colls":
						printf("\t\t\t\\(L\\) = <input type = \"text\" name = \"L\" size = \"3\" value = \"%s\"/>\n", isset($_POST["L"]) ? $_POST["L"] : "");
						printf("\t\t\t\\(S\\) = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
				}
				printf("\t\t</center>\n");
			}
		?>

		<?php
			if ($_POST["qty"] == "scatamp" or $_POST["qty"] == "dcs")
			{
				printf("\t\t<div class = \"text\">Set angular range:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(\\theta_{\mathrm{min}}\\) = <input type = \"text\" name = \"thmin\" size = \"5\" value = \"%s\"/>\n", $_POST["thmin"]);
				printf("\t\t\t\\(\\theta_{\mathrm{max}}\\) = <input type = \"text\" name = \"thmax\" size = \"5\" value = \"%s\"/>\n", $_POST["thmax"]);
				printf("\t\t\t\\(\\Delta\\theta\\) = <input type = \"text\" name = \"dth\" size = \"5\" value = \"%s\"/>\n", $_POST["dth"]);
				printf("\t\t</center>\n");
			}
		?>

		<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] == "ics" or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "colls"
				or $_POST["qty"] == "momtf" or $_POST["qty"] == "tcs")
			{
// 				// unit transform
// 				$Eunits = 1.;
// 				if (isset($_POST["Eunits"]))
// 				{
// 					switch ($_POST["Eunits"])
// 					{
// 						case "au":
// 							$Eunits = 0.5;
// 							break;
// 						case "eV":
// 							$Eunits = 13.605692;
// 							break;
// 						default:
// 							break;
// 					}
// 				}
				
				// default energies in Rydberg units
				$Emin_def = -1;
				$Emax_def = 0;
				$DE_def = 1;
				
				printf("\t\t<div class = \"text\" title = \"Set to '-1','0','1' to get all computed data. Otherwise you will get interpolated result. The interpolation is linear for most cases. Only for all integral cross sections at energies behind the ionization threshold the interpolation uses csplines.\">Set energy range:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(E_{\mathrm{min}}\\) = <input type = \"text\" name = \"Emin\" size = \"5\" value = \"%s\"/>\n", isset($_POST["Emin"]) ? $_POST["Emin"] : $Emin_def);
				printf("\t\t\t\\(E_{\mathrm{max}}\\) = <input type = \"text\" name = \"Emax\" size = \"5\" value = \"%s\"/>\n", isset($_POST["Emax"]) ? $_POST["Emax"] : $Emax_def);
				printf("\t\t\t\\(\\Delta E\\) = <input type = \"text\" name = \"dE\" size = \"5\" value = \"%s\"/>\n", isset($_POST["dE"]) ? $_POST["dE"] : $DE_def);
				printf("\t\t</center>\n");
			}
		?>
		
		<br/>
		
		<?php
			include "hexdbexe.inc";	// defines $hexdbexe
			include "hexdbdat.inc";	// defines $hexdbdat
			
			// units
			$strEunits = "Ry";
			if (isset($_POST["Eunits"]))
				$strEunits = $_POST["Eunits"];
			$strTunits = "a.u.";
			if (isset($_POST["Tunits"]))
				$strTunits = $_POST["Tunits"];
			
			if (isset($_POST["qty"]) and isset($_POST["view"]))
			{
				// prepare Hex-db command line
				$hexcmdline = $hexdbexe . " --database=" . $hexdbdat . " --" . $_POST["qty"];
				if (isset($_POST["ni"])) $hexcmdline = $hexcmdline . " --ni=" . $_POST["ni"];
				if (isset($_POST["li"])) $hexcmdline = $hexcmdline . " --li=" . $_POST["li"];
				if (isset($_POST["mi"])) $hexcmdline = $hexcmdline . " --mi=" . $_POST["mi"];
				if (isset($_POST["nf"])) $hexcmdline = $hexcmdline . " --nf=" . $_POST["nf"];
				if (isset($_POST["lf"])) $hexcmdline = $hexcmdline . " --lf=" . $_POST["lf"];
				if (isset($_POST["mf"])) $hexcmdline = $hexcmdline . " --mf=" . $_POST["mf"];
				if (isset($_POST["E"])) $hexcmdline = $hexcmdline . " --Ei=" . $_POST["E"];
				if (isset($_POST["L"])) $hexcmdline = $hexcmdline . " --L=" . $_POST["L"];
				if (isset($_POST["S"])) $hexcmdline = $hexcmdline . " --S=" . $_POST["S"];
				
				$hexcmdline = $hexcmdline . " --Eunits=" . $strEunits;
				$hexcmdline = $hexcmdline . " --Tunits=" . $strTunits;
				
				// compute standard input for Hex-db (energies or angles)
				if ($_POST["qty"] == "scatamp" or $_POST["qty"] == "dcs")
					$nums = range($_POST["thmin"], $_POST["thmax"], $_POST["dth"]);
				if ($_POST["qty"] == "ics" or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "colls"
					or $_POST["qty"] == "momtf" or $_POST["qty"] == "tcs")
					$nums = range($_POST["Emin"], $_POST["Emax"], $_POST["dE"]);

				// set PATH to include hex-db executable
//				putenv("PATH=" . "/home/jacob/Dokumenty/prog/Hex/hex-db/bin:" . $_ENV["PATH"]);
				
				// launch hex-db process
				$prochex = proc_open (
					$hexcmdline,
					array(array("pipe","r"), array("pipe","w"), array("pipe","a")),
					$pipes
				);
				
// 				stream_set_blocking($pipes[0], 0);
				
				// feed standard input to hex-db process
				foreach ($nums as $num)
					fwrite($pipes[0], $num . "\n");
				fclose($pipes[0]);
				
				// get standard output and close pipes
 				$hexoutput = stream_get_contents($pipes[1]);
				fclose($pipes[1]);
				fclose($pipes[2]);
				
				// terminate hex-db process
				$hex_return_value = proc_close($prochex);
				
				// store output data to hidden form element
				echo "<input type=\"hidden\" name=\"hexoutput\" value='$hexoutput' />";
			}
		?>

		<center>
			<input type = "submit" value = "View data" name = "view"/>
			&nbsp;&nbsp;&nbsp;
			<input type = "submit" value = "Download as TXT" name = "download" <?php
				if (!isset($_POST["qty"]) or !isset($_POST["view"]))
					echo "disabled=\"disabled\"";
			?>/>
		</center>

	</form>

	</td><td valign = "top" width = "50%">

	<div class = "sekce">Output:</div>
	<div class = "text">This section contains a graphical preview of selected data.</div>

	<?php
		if (isset($_POST["qty"]) and (isset($_POST["view"]) or isset($_POST["download"])))
		{
			// generate image
			$procgnuplot = proc_open(
				"/usr/bin/gnuplot",
				array(array("pipe","r"), array("pipe","w"), array("pipe","a")),
				$pipes2
			);
			
			// write to Gnuplot's standard input
			fwrite($pipes2[0], "set terminal png size 500,300\n");
			fwrite($pipes2[0], "unset key\n");
			fwrite($pipes2[0], "set xlabel \"Ei [" . $strEunits . "]\"\n");
			
			if ($_POST["qty"] == "colls")
				fwrite($pipes2[0], "set ylabel \"omega\"\n");
			else
				fwrite($pipes2[0], "set ylabel \"" . $_POST["qty"] . " [" . $strTunits . "]\"\n");
			
			if (isset($_POST["Tunits"]) and $_POST["Tunits"] == "cgs")
				fwrite($pipes2[0], "set format y '%g'\n");
				
// 			if ($nums[0] >= 1)
// 				fwrite($pipes2[0], "set logscale\n");
			
			if ($_POST["qty"] == "scatamp")
			{
				fwrite($pipes2[0], "set grid; plot [" . $nums[0] . ":" . end($nums) .  "] \"-\" using 1:2 with lines, \"\" using 1:3 with lines\n");
				fwrite($pipes2[0], $hexoutput);
				fwrite($pipes2[0], "e\n");
				fwrite($pipes2[0], $hexoutput);
				fwrite($pipes2[0], "e\n");
			}
			else
			{
				if ($_POST["Emin"] < 0 and $nums[0] < 0)
					fwrite($pipes2[0], "set grid; set logscale; plot \"-\" using 1:2 with lines\n");
				else
					fwrite($pipes2[0], "set grid; plot [" . $nums[0] . ":" . end($nums) .  "] \"-\" using 1:2 with lines\n");
				fwrite($pipes2[0], $hexoutput);
				fwrite($pipes2[0], "e\n");
			}
			fwrite($pipes2[0], "set terminal pop\n");
			
// 			stream_set_blocking($pipes2[0], 0);
			
			// close pipes
			fclose($pipes2[0]);
 			$gnuplot_out = stream_get_contents($pipes2[1]);
			fclose($pipes2[1]);
			fclose($pipes2[2]);
			
			// close the process
			$gnuplot_return_value = proc_close($procgnuplot);
			
			// display the plot
			echo "\t<div class = \"output\"><img src=\"data:image/png;base64," . base64_encode($gnuplot_out) . "\"/></div>\n";
			
			// write text data
// 			echo "\t<div class = \"output\"><pre>$hexoutput</pre></div>\n";
		}
	?>

	</td></tr><tr><td colspan="2" width = "100%">

		<div style="height: 1px; background-color: #880000; text-align: center">
			<span class = "sekce" style = "text-align: center; position: relative; top: -10px; background: white;">
				<a name = "avail-head"></a><a href = "#avail-head" id = "avail-link" onclick = "jsAvailClick()">&nbsp;&#x25BC; Available data &#x25BC;&nbsp;</a>
			</span>
		</div>
		<br/>
		<div class = "text" style = "display: none;" id = "avail-more-1">
			This section contains a graphical representation of data stored in the database
			at the moment. For every initial atomic state (vertical axis) the blue boxes
			show energy intervals (horizontal axis) covered by the data. In the (hopefully
			not so far) future all initial states will be covered completely with sufficient
			precision. Until that time one can use this chart as a simple measure
			of trusworthiness of the datasets. The darker the colour, the more partial waves
			have been included in the computation. A very simple way of how to verify that a
			particular chunk of energies has been computed with final precision is to compare
			"complete" and "extrapolated" cross section. If these two cross sections match,
			they ought to be reliable. The scattering amplitude may still not be converged,
			though, even in that case. The comparison is being served by default when
			the extrapolated cross section is requested: The resulting text file will
			contain both the "extrapolated" (\(\sigma_x\)) and the "complete" (\(\sigma_c\))
			cross section.
		</div>
		<div class = "output" style = "display: none;" id = "avail-more-2">
			<img src = "avail.png" alt = "avail.png"/>
		</div>
		<div class = "text" style = "display: none;" id = "avail-more-3">
			If your preview plot of a cross section contains a suspicious drop or rise,
			it may be a consequence of insufficient partial wave count. For the technical
			details on the computational settings that were used to produce the data see the
			<a href = "database.html">database</a> page.
		</div>

	</td></tr><tr><td colspan="2" width = "100%">

		<div style="height: 1px; background-color: #880000; text-align: center">
			<span class = "sekce" style = "text-align: center; position: relative; top: -10px; background: white;">
				<a name = "context-head"></a><a href = "#context-head" id = "context-link" onclick = "jsContextClick()">&nbsp;&#x25BC; Data in context &#x25BC;&nbsp;</a>
			</span>
		</div>
		<br/>
		<div class = "text" style = "display: none;" id = "context-more-1">
			Under construction... 
		</div>
		<div class = "output" style = "display: none;" id = "context-more-2">
			<img src = "empty.png" alt = "comparison"/>
		</div>

	</td></tr></table>

	</div> <!-- rám -->

	<div class = "pata">Jakub Benda &copy; 2013</div>

</body>

</html>
