#!/bin/csh -f
#  PlotAllStackModes.sh
#  RFHarmonicStacking
#
#  Created by Tolulope on 7/9/13.
#
#  Help understand how all stacks feed on each other ...
#  Currently, these stacks show epicentral and harmonic stacks.
#             I generalize this stacks for Azimuth harmonic / Azimuth unmodeeled harmonic

set rootDir = /Users/tmo22/Documents/OlugbojiPark2013/ParkOlugboji2012
set RFDebugDir = $rootDir/Code/DebugOutputs
set RFDataDir = $rootDir/SeismicData/SACbyStation/AllStationProcessed/RecFuncData
set Fig = $RFDebugDir/allStackTemplate.ps

#set MigStatus = "PpSmS Migration Depth"
#set MigStatus = "PpPmS Migration Depth"
set MigStatus = "Ps Migration Depth"
#set MigType = "MigrateReverb"
set MigType = "Migrate"
set freq = 1.2
set depth = 100
set width = 3
set height = 3
set bigHeight = 6
set offset = 0.1

set shwSynth = 0
set Cmpute = 1
set Mig = (No Yes Yes No)
set Azim = (No No Yes Yes)
set Epic = (No No No No)
set is = (2 3)

########################## BEGIN COMPUTE ALGORITHM >>>>>>>>>>>>>>>>>>>>>>>>>>>
if ($Cmpute) then
foreach i ($is)

	switch($Mig[$i])
		case [yY][eE][sS]:
			switch($Azim[$i])
				case [yY][eE][sS]:
					echo "$Mig[$i]" > DebugOptions
					echo "$Azim[$i]" >> DebugOptions
					echo "$Epic[$i]" >> DebugOptions
					echo "$depth" >> DebugOptions
					echo "$freq" >> DebugOptions
				    breaksw
				case [nN][oO]:
					echo "$Mig[$i]" > DebugOptions
					echo "$Azim[$i]" >> DebugOptions
					echo "$depth" >> DebugOptions
					echo "$freq" >> DebugOptions
				breaksw
			endsw
		breaksw
		case [nN][oO]:
				switch($Azim[$i])
					case [yY][eE][sS]:
						echo "$Mig[$i]" > DebugOptions
						echo "$Azim[$i]" >> DebugOptions
						echo "$Epic[$i]" >> DebugOptions
						echo "$freq" >> DebugOptions
					breaksw
					case [nN][oO]:
						echo "$Mig[$i]" > DebugOptions
						echo "$Azim[$i]" >> DebugOptions
						echo "$freq" >> DebugOptions
					breaksw
				endsw
		breaksw
	endsw

DebugOlugbojiRF.sh < DebugOptions >& /dev/null
echo $i

end
endif
######################## END COMPUTE ALGORITHM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

gmtset ANNOT_FONT_SIZE_PRIMARY 14p
set BOTLEFT = -JX$width"i"/$height"i"
set BOTRIGHT =  -JX$width"i"/$height"i"
			# Files in bottom columns ...
			set scale = 0.5
			set OutFilesHarm = $RFDebugDir/Debug
			set FileNoMig1 = $OutFilesHarm.Modelled.mean.xyz
			set FileNoMig2 = $OutFilesHarm.Modelled.devMin.xyz
			set FileNoMig3 = $OutFilesHarm.Modelled.devMax.xyz


			set OutFilesHarm = $RFDebugDir/Debug$MigType"Harmonic"
			set FileMig1 = $OutFilesHarm.Modelled.mean.xyz
			set FileMig2 = $OutFilesHarm.Modelled.devMin.xyz
			set FileMig3 = $OutFilesHarm.Modelled.devMax.xyz
			set Label = BazCmps


			# If  Migrate then show first second before migT
			set phase = `awk '(NR == 1){print $1}' $RFDebugDir/DebugtimeDelay.txt`
			set leftTime = ` (echo $phase) |awk '(NR == 1){print $1 - 2.0}' `
			set midTime = ` (echo $phase) |awk '(NR == 1){print $1 + 6.0}' `
			set rightTime = ` (echo $phase) |awk '(NR == 1){print $1 + 10.0}' `

			
			set BOXBOTLEFT = -R-2/10/-10/55
			set BOXBOTRIGHT = -R-2/10/-10/55
			#

set TOPLEFT = -JX$width"i"/$bigHeight"i"
set TOPRIGHT = -JX$width"i"/$bigHeight"i"

		# Files in top columns ...
		set scaleBot = 1.0
		set OutFilesHarm = $RFDebugDir/Debug
		set type = Epic
		set FileNoMigT1 = $OutFilesHarm.$type"_Rad.xyz"
		set FileNoMigT2 = $OutFilesHarm.$type"_Trans.xyz"

		# If Migrate ... Show Depth - Dist - RF
		set numFlds = `awk '{print NF}' $RFDebugDir/DebugMigrate.Epic_Trans.xyz | head -1`
		echo "*******************************" $numFlds
		if ( $numFlds == 4) then
			echo "Data Migrated With Depth ..."


			# Translate Time-Bounds to Depth Bounds ...
			set leftDepth = `awk '($2 >= -2.0){print $1}' $RFDebugDir/DebugMigrate.Epic_Rad.xyz | head -1`
			set rigthDepth = `awk '($2 >= 10.0){print $1}' $RFDebugDir/DebugMigrate.Epic_Rad.xyz | head -1`

			echo "Left " $leftDepth " Right    " $rigthDepth

			awk '(NF == 1){print $0}($2 >= -4 && $2 <= 15) {print $1 "  "  $3 "   " $4}' $RFDebugDir/DebugMigrate.Epic_Rad.xyz > $RFDebugDir/temp
			cp $RFDebugDir/temp $RFDebugDir/DebugMigrate.Epic_Rad_Dep.xyz
			#goto jumpout
		endif

		set FileNoMigS1 = $RFDataDir/"SynXMAS.Epic_Rad.xyz"
		set FileNoMigS2 = $RFDataDir/"SynXMAS.Epic_Trans.xyz"


		set OutFilesHarm = $RFDebugDir/Debug$MigType
		set FileMigT1 = $OutFilesHarm.$type"_Rad_Dep.xyz"
		set FileMigT2 = $OutFilesHarm.$type"_Trans.xyz"

		set FileMigS1 = $RFDataDir/"SynXMASMigrate.Epic_Rad.xyz"
		set FileMigS2 = $RFDataDir/"SynXMASMigrate.Epic_Trans.xyz"

		set Label = EpicentralDistance
		set BOXTOPLEFT = -R-2/10/25/105
		set BOXTOPRIGHT = -R-2/10/25/105
		set BOXTOPRIGHT = -R$leftDepth/$rigthDepth/25/105
		#



set Move = `echo $width " " $offset | awk '{print $1 + $2}'`

set MoveRight = -X$Move"i"
set MoveLeft = -X-$Move"i"

set Move = `echo $height " " $offset | awk '{print $1 + $2}'`

set MoveUp = -Y$Move"i"

#echo $TOPLEFT $TOPRIGHT $BOTLEFT $BOTRIGHT $MoveLeft $MoveRigth $MoveUp

# Draw BotLeft - Then Move Right
psbasemap $BOTLEFT $BOXBOTLEFT -Ba2f1g6:'Delay Time (sec)'::'.':/a70f10:'Baz Cmps ':/Swen -K -P >! $Fig

	pswiggle $FileNoMig3 $BOTLEFT $BOXBOTLEFT -Z$scale -M -G200/255/200  -V -P -O -K  >> $Fig
	pswiggle $FileNoMig2 $BOTLEFT $BOXBOTLEFT -Z$scale -M -G200/255/200 -N  -P -K -O  >> $Fig
	pswiggle $FileNoMig2 $BOTLEFT $BOXBOTLEFT -Z$scale -M -G0/0/255   -P -K -O  >> $Fig
	pswiggle $FileNoMig3 $BOTLEFT $BOXBOTLEFT -Z$scale -M -G255/0/0 -N  -P -K -O  >> $Fig

	pswiggle $FileNoMig1 $BOTLEFT $BOXBOTLEFT -Z$scale -W1p/0 -M  -P -K -O >>  $Fig

	echo -1.1 40 12 0 0 6 "Constant" | pstext $BOTLEFT -R-1/1/-10/50 -N -K -O  >> $Fig
	echo -1.13 30 12 0 0 6 "N-S" | pstext $BOTLEFT -R-1/1/-10/50 -N  -K -O  >> $Fig
	echo -1.13 20 12 0 0 6 "E-W" | pstext $BOTLEFT -R-1/1/-10/50 -N  -K -O  >> $Fig
	echo -1.13 10 12 0 0 6 "N-S" | pstext $BOTLEFT -R-1/1/-10/50 -N  -K -O  >> $Fig
	echo -1.13 0 12 0 0 6 "E-W" | pstext $BOTLEFT -R-1/1/-10/50 -N  -K -O  >> $Fig


# Draw BotRight - THen Move Up
	psbasemap $BOTRIGHT -R$leftTime/$rightTime/-10/55 -Ba2f1:'Delay Time (sec)'::'.':/a70:'Baz Cmps ':/Swen $MoveRight -O -K -P >> $Fig
	( echo $phase -10; echo $phase 55 ) | psxy $BOTRIGHT -R -O -K -P >> $Fig 
	( echo $midTime -10; echo $midTime 55 ) | psxy $BOTRIGHT -R -O -K -P >> $Fig

	pswiggle $FileMig3 $BOTRIGHT $BOXBOTRIGHT -Z$scale -M -G200/255/200     -V -O -K  >> $Fig
	pswiggle $FileMig2 $BOTRIGHT $BOXBOTRIGHT -Z$scale -M -G200/255/200 -N  -P -K -O  >> $Fig
	pswiggle $FileMig2 $BOTRIGHT $BOXBOTRIGHT -Z$scale -M -G0/0/255   -P -K -O  >> $Fig
	pswiggle $FileMig3 $BOTRIGHT $BOXBOTRIGHT -Z$scale -M -G255/0/0 -N  -P -K -O  >> $Fig

	pswiggle $FileMig1 $BOTRIGHT $BOXBOTRIGHT -Z$scale -W1p/0 -M  -P -K -O >>  $Fig


# Draw TopRight - Then Move Left
#psbasemap $TOPRIGHT -R -B:""::".":/:"":/swen $MoveUp -O -K -P >> $Fig


	pswiggle $FileMigT1 $TOPRIGHT $BOXTOPRIGHT $MoveUp -Z$scaleBot -M -G0/0/255 -Ba20f20:'Mig. Depth (km)'::'.':/a20f1:'Epicentral Distance (deg) ':/swEN   -V -O -K  >> $Fig
	( echo $depth  25; echo $depth 105 ) | psxy $TOPRIGHT -R -O -K -P >> $Fig
	pswiggle $FileMigT1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -M -G255/0/0 -N -B -P -K -O  >> $Fig

	pswiggle $FileMigT1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -W1p/0 -M  -P -K -O >>  $Fig

#echo 0 2.1 14 0 0 6 " $MigStatus $depth km " | pstext $TOPRIGHT -R-1/1/-2/2 -N -K -O  >> $Fig

if ( $shwSynth == 1) then
	pswiggle $FileMigS1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -M -G224/255/255  -P -K -O >>  $Fig
	pswiggle $FileMigS1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -M -G234/173/234 -N   -P -K -O >>  $Fig

	pswiggle $FileMigS1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -W1p/0 -M  -P -K -O >>  $Fig
endif

# Draw TopLeft - Then Rest in Peace.
#psbasemap $TOPLEFT -R -B:""::".":/:"":/swen $MoveLeft -O -K -P >> $Fig

if ( $shwSynth == 1) then
	pswiggle $FileNoMigS1 $TOPRIGHT $BOXTOPRIGHT  $MoveLeft -Z$scaleBot -M -G224/255/255  -P -K -O >>  $Fig
	pswiggle $FileNoMigS1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -M -G234/173/234 -N   -P -K -O >>  $Fig

	pswiggle $FileNoMigS1 $TOPRIGHT $BOXTOPRIGHT -Z$scaleBot -W1p/0 -M  -P -K -O >>  $Fig

	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT  -Z$scaleBot -M -G0/0/255 -Ba2f1g6:'Delay Time (sec)'::'.':/a20f5:'Epicentral Distance (deg)':/sWen   -V -O -K  >> $Fig
	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT -Z$scaleBot -M -G255/0/0 -N -B -P -K -O  >> $Fig

	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT -Z$scaleBot -W1p/0 -M  -P -K -O >>  $Fig

	echo 0 2.1 14 0 0 6 "XMAS | Freq $freq Hz | No Migration" | pstext $TOPLEFT -R-1/1/-2/2 -N -K -O  >> $Fig

else

	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT $MoveLeft -Z$scaleBot -M -G0/0/255 -Ba2f1g6:'Delay Time (sec)'::'.':/a20f5:'Epicentral Distance (deg)':/sWen   -V -O -K  >> $Fig
	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT -Z$scaleBot -M -G255/0/0 -N -B -P -K -O  >> $Fig

	pswiggle $FileNoMigT1 $TOPLEFT $BOXTOPLEFT -Z$scaleBot -W1p/0 -M  -P -K -O >>  $Fig

	echo 0 2.1 14 0 0 6 "XMAS | Freq $freq Hz | No Migration" | pstext $TOPLEFT -R-1/1/-2/2 -N -K -O  >> $Fig

endif





# View
open $Fig &

jumpout:


