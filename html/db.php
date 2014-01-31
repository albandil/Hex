<?php
	if (isset($_POST["download"]))
	{
		header("Content-Type: text/plain");
		header("Content-Disposition: attachment; filename=\"" . $_POST["qty"] . ".txt\"");
		echo $_POST["hexoutput"];
		exit();
	}

	// number of states to show in the table of available data (without ionization)
	$states = 9;
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
 
<head>
	<title>Hex scattering database</title>
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
	
	<link rel="icon" type="image/gif" href="hexe-small.gif" />
	
	<!-- Load stylesheet-->
	<style type="text/css">
		<!-- @import "style.css"; -->
	</style>

	<!-- Copy value of $states [from PHP] to Nstates [JavaScript] -->
	<script type = "text/javascript" language = "javascript">
		<?php echo "Nstates = $states\n"; ?>
	</script>
	
	<!-- Load external scripts -->
	<script type = "text/javascript" language = "javascript" src = "db-scripts.js">
	</script>

	<!-- Setup MathJax -->
	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
			extensions: ["tex2jax.js"],
			jax: ["input/TeX","output/HTML-CSS"],
// 			tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
		});
	</script>

	<!-- Load MathJax -->
	<!--http://www.mathjax.org/mathjax/MathJax.js-->
	<script src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
</head>
 
<body onload = "jsDataSet(); jsDataAngular();">

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

	<!-- body -->

	<table width = "100%" style = "table-layout:fixed;"><tr><td valign = "top" width = "50%">

	<div class = "sekce">Input:</div>

	<!-- scattering variable -->
	<form name = "data" action = "db.php" method = "post" style = "margin: 10px">
		<div class = "text">Choose what to compute:</div>
		<center>
			<select name = "qty" title = "scattering quantity to compute" onchange = "this.form.submit()">
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

		<!-- units -->
		<div class = "text">Choose units:</div>
		<center>
			<span class = "text">Energy: <select name = "Eunits" title = "energy units for energy input">
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
			</select> ,  Output: <select name = "Tunits" title = "length units for (dimensioned) output">
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

		<!-- initial atomic state -->
		<div class = "text">Set initial atomic state(s):</div>
		<center>
			\(n_i\) = <input type = "text" title = "initial principal quantum number" name = "ni" size = "3" value = "<?php echo (isset($_POST["ni"]) ? $_POST["ni"] : 1); ?>" required = "required"/>
			\(l_i\) = <input type = "text" title = "initial orbital quantum number" name = "li" size = "3" value = "<?php echo (isset($_POST["li"]) ? $_POST["li"] : 0); ?>" required = "required"/>
			\(m_i\) = <input type = "text" title = "initial magnetic quantum number" name = "mi" size = "3" value = "<?php echo (isset($_POST["mi"]) ? $_POST["mi"] : 0); ?>" required = "required"/>
		</center>
		
		<!-- final tomic state -->
<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] != "tcs")
			{
				printf("\t\t<div class = \"text\">Set final atomic state(s):</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(n_f\\) = <input type = \"text\" title = \"final principal quantum number\" name = \"nf\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["nf"]) ? $_POST["nf"] : "1");
				printf("\t\t\t\\(l_f\\) = <input type = \"text\" title = \"final orbital quantum number\" name = \"lf\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["lf"]) ? $_POST["lf"] : "0");
				printf("\t\t\t\\(m_f\\) = <input type = \"text\" title = \"final magnetic quantum number\" name = \"mf\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["mf"]) ? $_POST["mf"] : "0");
				printf("\t\t</center>\n");
			}
?>
		
		<!-- total quantum numbers -->
<?php
			if (isset($_POST["qty"]) and $_POST["qty"] == "scatamp") //~ scatamp
			{
				printf("\t\t<div class = \"text\">Set total quantum numbers:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(E\\) = <input type = \"text\" title = \"impact energy of the incoming electron\" name = \"E\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
				printf("\t\t\t\\(S\\) = <input type = \"text\" title = \"total spin of the two electrons\" name = \"S\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
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
						printf("\t\t\t\\(E\\) = <input type = \"text\" title = \"impact energy of the incoming electron\" name = \"E\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
						printf("\t\t\t\\(S\\) = <input type = \"text\" title = \"total spin of the two electrons\" name = \"S\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "momtf":
						printf("\t\t\t\\(S\\) = <input type = \"text\" title = \"total spin of the two electrons\" name = \"S\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "ics":
					case "colls":
						printf("\t\t\t\\(L\\) = <input type = \"text\" title = \"total orbital momentum of the two electrons\" name = \"L\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["L"]) ? $_POST["L"] : "");
						printf("\t\t\t\\(S\\) = <input type = \"text\" title = \"total spin of the two electrons\" name = \"S\" size = \"3\" value = \"%s\" required = \"required\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
				}
				printf("\t\t</center>\n");
			}
?>

		<!-- scattering angles -->
<?php
			if ($_POST["qty"] == "scatamp" or $_POST["qty"] == "dcs")
			{
				printf("\t\t<div class = \"text\">Set angular range:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(\\theta_{\mathrm{min}}\\) = <input type = \"text\" title = \"smallest scattering angle\" name = \"thmin\" size = \"5\" value = \"%s\" required = \"required\"/>\n", $_POST["thmin"]);
				printf("\t\t\t\\(\\theta_{\mathrm{max}}\\) = <input type = \"text\" title = \"largest scattering angle\" name = \"thmax\" size = \"5\" value = \"%s\" required = \"required\"/>\n", $_POST["thmax"]);
				printf("\t\t\t\\(\\Delta\\theta\\) = <input type = \"text\" title = \"spacing between the scattering angles\" name = \"dth\" size = \"5\" value = \"%s\" required = \"required\"/>\n", $_POST["dth"]);
				printf("\t\t</center>\n");
			}
?>

		<!-- impact energies -->
<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] == "ics" or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "colls"
				or $_POST["qty"] == "momtf" or $_POST["qty"] == "tcs")
			{
				// get checkbox status
				if (isset($_POST["Eall"]) or !isset($_POST["qty"]))
				{
					$editstatus = " disabled = \"disabled\"";
					$chckstatus = " checked = \"checked\"";
				}
				else
				{
					$editstatus = "";
					$chckstatus = "";
				}
				
				printf("\t\t<div class = \"text\" title = \"Set to '-1','0','1' to get all computed data. Otherwise you will get interpolated result. The interpolation is linear for most cases. Only for all integral cross sections at energies behind the ionization threshold the interpolation uses csplines.\">Set energy range:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t\\(E_{\mathrm{min}}\\) = <input type = \"text\" title = \"lowest impact energy\" id = \"iEmin\" name = \"Emin\" size = \"5\" value = \"%s\" $editstatus required = \"required\"/>\n", isset($_POST["Emin"]) ? $_POST["Emin"] : "");
				printf("\t\t\t\\(E_{\mathrm{max}}\\) = <input type = \"text\" title = \"highest impact energy\" id = \"iEmax\" name = \"Emax\" size = \"5\" value = \"%s\" $editstatus required = \"required\"/>\n", isset($_POST["Emax"]) ? $_POST["Emax"] : "");
				printf("\t\t\t\\(\\Delta E\\) = <input type = \"text\" title = \"impact energy spacing\" id = \"idE\" name = \"dE\" size = \"5\" value = \"%s\" $editstatus required = \"required\"/>\n", isset($_POST["dE"]) ? $_POST["dE"] : "");
				printf("\t\t</center>\n");
				
				
				printf("\t\t<div class = \"text\">or</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\t<input type = \"checkbox\" name = \"Eall\" onclick = \"jsEallClick()\" value = \"1\"$chckstatus> retrieve all available energies.\n");
				printf("\t\t</center>\n");
			}
?>
		
		<br/>
		
		<!-- hidden element containing the output from hex-db -->
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
				
				if (isset($_POST["Eall"]))
				{
					$_POST["Emin"] = -1;
					$_POST["Emax"] =  0;
					$_POST["dE"]    = 1;
				}
				
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

		<!-- view/download buttons -->
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
	<div class = "text">This section contains a graphical preview of the selected data.</div>

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
// 			fwrite($pipes2[0], "set terminal svg mouse jsdir \"http://gnuplot.sourceforge.net/demo_svg_4.6/\" size 500,300\n"); // SVG
			fwrite($pipes2[0], "set terminal png size 500,300\n"); // PNG
			fwrite($pipes2[0], "unset key\n");
			fwrite($pipes2[0], "set xlabel \"Ei [" . $strEunits . "]\"\n");
			
			if ($_POST["qty"] == "colls")
				fwrite($pipes2[0], "set ylabel \"omega\"\n");
			else
				fwrite($pipes2[0], "set ylabel \"" . $_POST["qty"] . " [" . $strTunits . "]\"\n");
			
			if (isset($_POST["Tunits"]) and $_POST["Tunits"] == "cgs")
				fwrite($pipes2[0], "set format y '%g'\n");
				
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
			
			// display the SVG plot
// 			echo $gnuplot_out;
			
			// display the PNG plot
			echo "\t<div class = \"output\"><img src=\"data:image/png;base64," . base64_encode($gnuplot_out) . "\"/></div>\n";
			
			// display input data for debugging
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
		<div class = "text" style = "display:none;" id = "avail-more-1">
			The simple table below ilustrates the current state of the contents
			of the database. The rows of the table are different initial atomic
			states (before the collision), the columns are different final
			states. Select particular angular momenta from the drop-down menus,
			pick a cell that corresponds to the principal quantum number and
			the available data will be shown. The colours in the table represent
			subjective rating of the data for the specific transition. White colour
			stands for "no data at all". Dark colour means lots of angular momentum
			transitions, for lots of energies. If you need some specific data
			that are not present, do not hesitate to contact the author.
		</div>
		
		<center><table class = "availdata" id = "avail-more-2" style = "display:none;">
			<colgroup>
				<col/>
				<?php $W = 55/($states + 2); for ($i = 0; $i <= $states+1; $i++) echo "<col width=\"$W%\">"; ?>
				<col/>
			</colgroup>
			<tr>
				<td rowspan = "2" colspan = "2" style = "border-width:0;"></td>
				<td colspan = "<?php echo ($states+1); ?>" align = "center" style = "border-width:0;">
					\(n_f\)
					<select id = "alf" title = "lf" onchange = "jsDataAngular()">
						<option value = "s" selected = "selected">s</option>
						<option value = "p">p</option>
						<option value = "d">d</option>
						<option value = "f">f</option>
						<option value = "g">g</option>
						<option value = "h">h</option>
						<option value = "i">i</option>
						<option value = "j">j</option>
					</select>
				</td>
				<td style = "border-width:0;"></td>
			</tr>
			<tr>
				<?php for ($i = 1; $i <= $states; $i++) echo "<td id = \"head-f-$i\" bgcolor = \"#880000\" style = \"color: white;\">$i</td>"; ?>
				<td bgcolor = "#880000" style = "color: white;">ion.</td>
			</tr>
			<tr>
				<td rowspan = "<?php echo $states; ?>" align = "center" style = "border-width:0;">
					\(n_i\)
					<br/>
					<select id = "ali" title = "li" onchange = "jsDataAngular()">
						<option value = "s" selected = "selected">s</option>
						<option value = "p">p</option>
						<option value = "d">d</option>
						<option value = "f">f</option>
						<option value = "g">g</option>
						<option value = "h">h</option>
						<option value = "i">i</option>
						<option value = "j">j</option>
					</select>
				</td>
				<td id = "head-i-1" bgcolor = "#880000" style = "color: white;">1</td>
				<?php for ($i = 1; $i <= $states; $i++) echo "<td id = \"dat-1-$i\"></td>"; echo "\n"; ?>
				<td id = "dat-1-0"></td>
				<td rowspan = "<?php echo $states; ?>" width = "40%" id = "datadescr" valign = "top">
					<!-- notes, to be written by JS -->
				</td>
			</tr>
<?php
			for ($i = 2; $i <= $states; $i++)
			{
				echo "\t\t\t<tr>\n";
				echo "\t\t\t\t<td id = \"head-i-$i\" bgcolor = \"#880000\" style = \"color: white;\">$i</td>";
				for ($j = 1; $j <= $states; $j++)
					echo "<td id = \"dat-$i-$j\"></td>";
				echo "<td id = \"dat-$i-0\"></td>";
				echo "\n\t\t\t</tr>\n";
			}
?>
		</table></center>
	
	</td></tr></table>

	</div> <!-- rám -->

	<div class = "pata"><a href = "mailto:jakub.benda@seznam.cz?subject=Hex web">Jakub Benda</a> &copy; 2014</div>

</body>

</html>
