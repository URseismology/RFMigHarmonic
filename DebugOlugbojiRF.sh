#!/bin/csh -f
# Updated on Nov. 14, 2012
# script to Debug and Iteratively develop the Harmonic stack mode: -M option
#    Look at Code Blog to follow thoughts....
# Currently built to implement methods in JJParks 2010 analysis and legacy code rfmig_cboot.f
#
# From legacy code to version control in git. Repos. And all that.

# Output Files returned on execution of script
#                     1. Azimuth_Rad.xyz     Azimuth_Trans.xyz
#                     2. OutputFName.log      : Stack & Bin Statistics
#                     3. OutputFNameEvent.txt : Event Statistics
#											  : Azimuth, EvLong,Evlat,Mag,Dist.

# Never Push this Code to the repo.
# Either I always put out two different
set onMac = 1

if ($onMac) then

    set rootDir = /Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012
    set Distrib = $rootDir/Code/RFMigHarmonic
    set RecordDir = $rootDir/SeismicData/processedSeismograms/StatFiles
    set OutFiles = $rootDir/Code/DebugOutputs
    set Exec = $rootDir/Code/bin
	set SACDir = $rootDir/SeismicData/SACbyStation/AllStationProcessed/SACFiles
	set RFDataDir = $rootDir/SeismicData/SACbyStation/AllStationProcessed/RecFuncData
	set LocalMake = $Distrib/MakefileOnMAC


else
    # Code here is used to Run on LAB Machines. This allows me to have a dual debug
    # Interface.
	set rootDir = /Users/tmo22/Documents/ParkOlugboji2012/
	set Distrib = /Users/tmo22/Documents/ParkOlugboji2012/CodeRepo/RFMigHarmonic
	set RecordDir = /Users/tmo22/Documents/ParkOlugboji2012/SeismicData/processedSeismograms/StatFiles
	set OutFiles = /Users/tmo22/Documents/ParkOlugboji2012/CodeRepo/DebugOutputs
	set SACDir = $rootDir/SeismicData/SACbyStation/AllStationProcessed/SACFiles
	set RFDataDir = $rootDir/SeismicData/SACbyStation/AllStationProcessed/RecFuncData
	set Exec = /Users/tmo22/Documents/ParkOlugboji2012/CodeRepo/bin

	set LocalMake = $Distrib/Makefile

endif

set migrateFile = $rootDir/Code/ShellDataPrep/XMAS_IsoHKmodels.txt

set isDisplay = 1					# Set If Files Are to Be Displayed
#set stationManifest = XMAS.in
set stationManifest = XMASMAC.in
#set stationManifest = DataDebug.in

##################################### Wire forward modelling into this shell script.
set runFwdMdels = 0

if ($runFwdMdels) then
	$rootDir/Code/ShellDataPrep/forwardmodelHarmonic.sh
endif

set runFwdMdels = 0

######## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#goto showResults

# I use DataDebug.in as a record list for testing code functions.
# I  have two RF modules - RF data and RFH-K modelling ...
#

set isRForHK = 1   # If RF then Compile RecFunc, Otherwise Compile HK Stacking Code

if ($isRForHK == 1 ) then
	make help -f $LocalMake
	make RecFunc -f $LocalMake
	cp $Distrib/RecFunc $Exec/RecFuncHarmonicJack
else
	make help -f $LocalMake
	make RFVelStck -f $LocalMake
	
	### Independent tests for now... Working on Sequential HkVk Stacks
	set makeStatus = $status

	if ( $makeStatus == 2 ) then
		echo "Error in Make file revisit its execution"
	else
		echo  "Compilation Success... Running Basics"
		set srchParamsFile = $rootDir/Code/GRDSrchParams/KIP_GRDsrch.txt
		$Distrib/RFVelStck -L$RecordDir/$stationManifest -O$OutFiles/RFVelDebug -S$srchParamsFile -R0.75 -I5
	endif

	goto finish
	
endif

set makeStatus = $status

if ( $makeStatus == 2 ) then
	echo "Error in Make file revisit its execution"
endif

goto finish
####################################  RF Parameters Set Here >>>>>>>>>>>>>>>>>>>>>>>>>>>
set leftAzim = 0
set rightAzim = 360
set epicStep = 2
set rotate = 0.65
set migrationFlag =	-M0/				# Direct Phase Reverberation
#set migrationFlag =	-m0/				# Reverberated Phase - PpSms
#set migrationFlag =	-m1/				# Reverberated Phase - PpPms
# END RF Parameter Set Here >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cd $SACDir							 # Check this...
if ( $makeStatus == 0 ) then
	cp $Distrib/RecFunc $Exec/RecFuncHarmonic

	echo "Migrate? Yes or No"
	set migrateResponse = $<

	switch ($migrateResponse)
		case [yY][eE][sS]:
			echo "Running Migration.. .Please Provide More Responses ..."
			echo "Azimuth-Epicetral[AE] or Harmonic[H]?. Yes for [AE] or No for [H]"
			set stackResponse = $<
				switch ($stackResponse)
					case [yY][eE][sS]:
						echo "Azimuth or Epicentral? Yes for Azimuth Stack. If No, then Epicentral Stack"
						set isAzimuth = $<

						switch ($isAzimuth)
							case [yY][eE][sS]:
								# Azimuth Stack
								echo "Running Azimuth Stack - Enter Migration Depth in Km"
								set migDepth = $<

								echo "Enter Freq. in Hz..."
								set freq = $<
								
								#####
								$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
									-F$freq -T40 -O$OutFiles/Debug  -A0/350/5 \
									-B3 -H1 -V -S0/50 -R$rotate \
									-I$migrateFile $migrationFlag$migDepth
								
								##### RF for Synthetics 
							if ($runFwdMdels) then	

								$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth \
								-F$freq -T40 -O$RFDataDir/SynXMAS  -A0/350/5 \
								-B3 -H1 -V -S0/50 -R$rotate \
								-I$migrateFile $migrationFlag$migDepth
							endif

								breaksw
							case [nN][oO]:
								# Epicentral Stack
								echo "Running Epicentral Stack - Enter Migration Depth in Km"
								set migDepth = $<

								echo "Enter Freq. in Hz..."
								set freq = $<

								$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
								-F$freq -T40 -O$OutFiles/Debug \
								-E$leftAzim/$rightAzim/30/100/$epicStep -B1 -H1 -V \
								-S0/50 -R$rotate \
								-I$migrateFile $migrationFlag$migDepth

							##### RF for Synthetics 
							if ($runFwdMdels) then	
								$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth\
								-F$freq -T40 -O$RFDataDir/SynXMAS \
								-E$leftAzim/$rightAzim/30/100/$epicStep -B1 -H1 -V \
								-S0/50 -R$rotate \
								-I$migrateFile $migrationFlag$migDepth
							endif

								breaksw
							default:
								echo "Can't determine if Azimuth or Epicentral..."
								breaksw
						endsw
						breaksw

					case [nN][oO]:
						echo "Running Harmonic Stack - Enter Migration Depth in Km"
						set migDepth = $<

						echo "Enter Freq. in Hz..."
						set freq = $<

						# Harmonic Stack
						$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
							-F$freq -T40 -O$OutFiles/Debug -H1 -V -S1/100 -R0.75 \
							-I$migrateFile $migrationFlag$migDepth

						##### RF for Synthetics 
						if ($runFwdMdels) then	
							$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth\
							-F$freq -T40 -O$RFDataDir/SynXMAS -H1 -V -S1/100 -R0.75 \
							-I$migrateFile $migrationFlag$migDepth
						endif

						breaksw
					default:
						echo "Can't determine what stack mode is set ..."
				endsw
			breaksw

		case [nN][oO]:
			echo "Not Running Migration ... Please Provide More Responses ..."

			echo "Azimuth Stack or Harmonic Stack. Yes or No? if No, then HarmonicStack"
			set stackResponse = $<
			switch ($stackResponse)
				case [yY][eE][sS]:
					echo "Running Azim or/and Harmonic Stack ..."
					echo "Answer. Yes for Azimuth Stack. If No, then Epicentral Stack"
					set isAzimuth = $<

					echo "Enter Freq. in Hz..."
					set freq = $<

					switch ($isAzimuth)
						case [yY][eE][sS]:
							# Azimuth Stack
							$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
								-F$freq -T40 -O$OutFiles/Debug  -B3 -A0/350/5 \
								-H1 -V -S0/50 -R$rotate

							##### RF for Synthetics 
							if ($runFwdMdels) then	
								$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth \
								-F$freq -T40 -O$RFDataDir/SynXMAS  -B3 -A0/350/5 \
								-H1 -V -S0/50 -R$rotate
							endif

							breaksw
						case [nN][oO]:
							# Epicentral Stack
							$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
								-F$freq -T40 -O$OutFiles/Debug  -B1 \
								-E$leftAzim/$rightAzim/30/100/$epicStep -H1 -V -S0/50 -R$rotate

							##### RF for Synthetics 
							if ($runFwdMdels) then	
								$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth \
								-F$freq -T40 -O$RFDataDir/SynXMAS  -B1 \
								-E$leftAzim/$rightAzim/30/100/$epicStep -H1 -V -S0/50 -R$rotate
							endif

							breaksw
						default:
							echo "Can't determine if Azimuth or Epicentral..."
							breaksw
				
					endsw

					breaksw
				case [nN][oO]:
					echo "Running Harmonic Stack"

					echo "Enter Freq. in Hz..."
					set freq = $<
					# Harmonic Stack
					$Exec/RecFuncHarmonic -L$RecordDir/$stationManifest \
						-F$freq -T40 -O$OutFiles/Debug  -B3 -H1 -V -S1/50 -R0.65

					##### RF for Synthetics 
					if ($runFwdMdels) then	
						$Exec/RecFuncHarmonic -L$SACDir/in_recpick_synth \
						-F$freq -T40 -O$RFDataDir/SynXMAS  -B3 -H1 -V -S1/50 -R0.65
					endif

					breaksw
				default:
					echo "Can't determine what stack mode is set ..."
			endsw
			breaksw
		default:
			echo "Can't determine what migration option is set."
	endsw

	mv *.xyz $OutFiles/

#-I$rootDir/Code/ShellDataPrep/XMAS_AniMod.txt -M1/50
endif


goto finish
showResults:

if ( $isDisplay == 1) then
set RFDataDir = /Users/tmo22/Documents/ParkOlugboji2012/Code/bin
set FigOut = /Users/tmo22/Documents/ParkOlugboji2012/Images

echo "EDIT PLOTTING PARAMETERS TO SUIT YOUR DATA"

	if (0) then
		set FRAME = -JX3.35i/6i
		set FRAME_BIG = -JX3i/7i
		set BOX = -R-4/8/-10/55

		set scale = 0.3
		set OutFilesHarm = $OutFiles/Debug"MigrateHarmonic"
		set File1 = $OutFilesHarm.Modelled.mean.xyz
		set File2 = $OutFilesHarm.UnModelled.mean.xyz
		set Label = BazCmps
		set Station = Debug$stationManifest
	else
		set FRAME = -JX3.35i/6i
		set FRAME_BIG = -JX3i/7i
		#set BOX = -R-4/8/20/110

		set BOX = -R-4/8/0/380		# Azimuth
		set type = Azim 
		set scale = 0.8
		set File1 = $OutFiles/DebugMigrate.$type"_Rad.xyz"
		set File2 = $OutFiles/DebugMigrate.$type"_Trans.xyz"
		set Label = BazCmps
		set Station = Debug$stationManifest"Migrate"
	endif


	set Fig = $OutFiles/$Station.ps

	pswiggle $File1 $FRAME $BOX -Z$scale -M -G0/0/255 -Ba2f1g8:'Delay Time (sec)'::'.Modelled':/a20f1:'Baz Cmps ':/SWen -Y4.0  -V -P -K  >! $Fig
	pswiggle $File1 $FRAME $BOX -Z$scale -M -G255/0/0 -N -B -P -K -O  >> $Fig

	pswiggle $File1 $FRAME $BOX -Z$scale -W1p/0 -M  -P -K -O >>  $Fig
	#echo 0 2.7 14 0 5 6 $File1 scale $scale | pstext $FRAME -R-1/1/-2/2 -N -O -K >> $Fig


	pswiggle $File2 $FRAME $BOX -Z$scale -M -G0/0/255 -Ba2f1:'Delay Time (sec)'::'.UnModelled':/a20f10:'':/SwEn -V -X3.35i -P -O -K  >> $Fig
	pswiggle $File2 $FRAME $BOX -Z$scale -M -G255/0/0 -N -B -P -K -O  >> $Fig

	pswiggle $File2 $FRAME $BOX -Z$scale -W1p/0 -M  -P -K -O >>  $Fig
	echo 0 2.7 14 0 0 6 $Station scale $scale | pstext $FRAME -R-1/1/-2/2 -N -O  >> $Fig

	open $Fig &

endif

finish:
