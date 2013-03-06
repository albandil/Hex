<?php

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *                       / /   / /    __    \ \  / /                         *
 *                      / /__ / /   / _ \    \ \/ /                          *
 *                     /  ___  /   | |/_/    / /\ \                          *
 *                    / /   / /    \_\      / /  \ \                         *
 *                                                                           *
 *                         Jakub Benda (c) 2013                              *
 *                     Charles University in Prague                          *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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

		<div class = "text">Set initial atomic state(s):</div>
		<center>
			ni = <input type = "number" name = "ni" size = "3" value = "<?php echo (isset($_POST["ni"]) ? $_POST["ni"] : 1); ?>"/>
			li = <input type = "number" name = "li" size = "3" value = "<?php echo (isset($_POST["li"]) ? $_POST["li"] : 0); ?>"/>
			mi = <input type = "number" name = "mi" size = "3" value = "<?php echo (isset($_POST["mi"]) ? $_POST["mi"] : 0); ?>"/>
		</center>
		
		<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] != "tcs")
			{
				printf("\t\t<div class = \"text\">Set final atomic state(s):</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\tnf = <input type = \"text\" name = \"nf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["nf"]) ? $_POST["nf"] : "1");
				printf("\t\t\tlf = <input type = \"text\" name = \"lf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["lf"]) ? $_POST["lf"] : "0");
				printf("\t\t\tmf = <input type = \"text\" name = \"mf\" size = \"3\" value = \"%s\"/>\n", isset($_POST["mf"]) ? $_POST["mf"] : "0");
				printf("\t\t</center>\n");
			}
		?>
		
		<?php
			if (isset($_POST["qty"]) and $_POST["qty"] == "scatamp") //~ scatamp
			{
				printf("\t\t<div class = \"text\">Set global quantum numbers:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\tE = <input type = \"text\" name = \"E\" size = \"3\" value = \"%s\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
				printf("\t\t\tS = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
				printf("\t\t</center>\n");
			}
			else if (!isset($_POST["qty"]) or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "tcs")
			{
				// do nothing
			}
			else
			{
				printf("\t\t<div class = \"text\">Set global quantum numbers:</div>\n");
				printf("\t\t<center>\n");
				switch ($_POST["qty"])
				{
					case "scatamp":
					case "dcs":
						printf("\t\t\tE = <input type = \"text\" name = \"E\" size = \"3\" value = \"%s\"/>\n", isset($_POST["E"]) ? $_POST["E"] : "");
						printf("\t\t\tS = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "momtf":
						printf("\t\t\tS = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
						break;
					case "ics":
					case "colls":
						printf("\t\t\tL = <input type = \"text\" name = \"L\" size = \"3\" value = \"%s\"/>\n", isset($_POST["L"]) ? $_POST["L"] : "");
						printf("\t\t\tS = <input type = \"text\" name = \"S\" size = \"3\" value = \"%s\"/>\n", isset($_POST["S"]) ? $_POST["S"] : "");
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
				printf("\t\t\tθmin = <input type = \"text\" name = \"thmin\" size = \"5\" value = \"%s\"/>\n", $_POST["thmin"]);
				printf("\t\t\tθmax = <input type = \"text\" name = \"thmax\" size = \"5\" value = \"%s\"/>\n", $_POST["thmax"]);
				printf("\t\t\tΔθ = <input type = \"text\" name = \"dth\" size = \"5\" value = \"%s\"/>\n", $_POST["dth"]);
				printf("\t\t</center>\n");
			}
		?>

		<?php
			if (!isset($_POST["qty"]) or $_POST["qty"] == "ics" or $_POST["qty"] == "ccs" or $_POST["qty"] == "xcs" or $_POST["qty"] == "colls"
				or $_POST["qty"] == "momtf" or $_POST["qty"] == "tcs")
			{
				printf("\t\t<div class = \"text\">Set energy range:</div>\n");
				printf("\t\t<center>\n");
				printf("\t\t\tEmin = <input type = \"text\" name = \"Emin\" size = \"5\" value = \"%s\"/>\n", isset($_POST["Emin"]) ? $_POST["Emin"] : 0.65);
				printf("\t\t\tEmax = <input type = \"text\" name = \"Emax\" size = \"5\" value = \"%s\"/>\n", isset($_POST["Emax"]) ? $_POST["Emax"] : 0.85);
				printf("\t\t\tΔE = <input type = \"text\" name = \"dE\" size = \"5\" value = \"%s\"/>\n", isset($_POST["dE"]) ? $_POST["dE"] : 0.001);
				printf("\t\t</center>\n");
			}
		?>
		
		<br/>
		
		<?php
			include "hexdbexe.inc";	// defines $hexdbexe
			include "hexdbdat.inc";	// defines $hexdbdat
			
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

	<div class = sekce>Output:</div>

	<?php
		if (isset($_POST["qty"]) and (isset($_POST["view"]) or isset($_POST["download"])))
		{
			// generate image
			$procgnuplot = proc_open(
				"/usr/bin/gnuplot",
				array(array("pipe","r"), array("pipe","w"), array("pipe","a")),
				$pipes2
			);
			
			// write to Gnuplot standard input
			fwrite($pipes2[0], "set terminal png size 500,300\n");
			fwrite($pipes2[0], "unset key\n");
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
		<center>
			<div class = "sekce">Available data:</div>
			<div class = "output"><img src = "avail.png" alt = "avail.png"/></div>
		</center>
	</td></tr></table>

	</div> <!-- rám -->

	<div class = "pata">Jakub Benda ⓒ 2013</div>

</body>

</html>
