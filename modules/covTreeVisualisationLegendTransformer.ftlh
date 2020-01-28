<#setting number_format="####0.###">
<#assign legendVisDoc = transformerInput.treeVisualisationLegend>
<#assign svgWidthPx = legendVisDoc.legendPxWidth>
<#assign svgHeightPx = legendVisDoc.legendPxHeight>
<#assign legendFontSize = legendVisDoc.legendFontSize>
	
<svg width="${svgWidthPx}" height="${svgHeightPx}"
     xmlns="http://www.w3.org/2000/svg">
  <style>${getResourceAsString('covTreeVisualisationTransformer.css')?no_esc}</style>
  <g class="branchLengthLegendLine">
    <line
	  	x1="${legendVisDoc.branchLengthLegendCentreLine.x1}" 
	  	y1="${legendVisDoc.branchLengthLegendCentreLine.y1}" 
	  	x2="${legendVisDoc.branchLengthLegendCentreLine.x2}" 
	  	y2="${legendVisDoc.branchLengthLegendCentreLine.y2}"/>
    <#list legendVisDoc.branchLengthLegendTicks as branchLengthLegendTick>
    <line
	  	x1="${branchLengthLegendTick.x1}" 
	  	y1="${branchLengthLegendTick.y1}" 
	  	x2="${branchLengthLegendTick.x2}" 
	  	y2="${branchLengthLegendTick.y2}"/>
    </#list>
  </g>    
  <g class="legendText" font-size="${legendFontSize}">
  	<#list legendVisDoc.branchLengthLegendValues as branchLengthLegendValue>
    <text
	  	x="${branchLengthLegendValue.x}"
	  	y="${branchLengthLegendValue.y}"
	  	width="${branchLengthLegendValue.width}" 
	  	height="${branchLengthLegendValue.height}"
	  	dy=".35em">${branchLengthLegendValue.text}</text>
  	</#list>
    <text
	  	x="${legendVisDoc.branchLengthLegendUnits.x}"
	  	y="${legendVisDoc.branchLengthLegendUnits.y}"
	  	width="${legendVisDoc.branchLengthLegendUnits.width}" 
	  	height="${legendVisDoc.branchLengthLegendUnits.height}"
	  	dy=".35em">${legendVisDoc.branchLengthLegendUnits.text}</text>
  </g>
</svg>
