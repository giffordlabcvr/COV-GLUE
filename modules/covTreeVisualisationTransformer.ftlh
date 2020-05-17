<#setting number_format="####0.###">
<#assign treeVisDoc = transformerInput.treeVisualisation>
<#assign svgWidthPx = treeVisDoc.pxWidth>
<#assign svgHeightPx = treeVisDoc.pxHeight>
<#assign verticalLeafSpacePx = treeVisDoc.verticalLeafSpacePx>
<#assign leafTextGapPx = treeVisDoc.leafTextGapPx>
<#assign leafTextFontSize = treeVisDoc.leafTextFontSize>
	
<#assign leafNodeRadiusPx = verticalLeafSpacePx * 0.6>
	
<svg version="1.1" width="${svgWidthPx}" height="${svgHeightPx}"
     xmlns="http://www.w3.org/2000/svg">
  <style>${getResourceAsString('covTreeVisualisationTransformer.css')?no_esc}</style>

  <linearGradient id="collapsedGrad" x1="0%" y1="0%" x2="100%" y2="0%">
    <stop offset="0%" style="stop-color:DarkGrey;stop-opacity:1" />
    <stop offset="60%" style="stop-color:DarkGrey;stop-opacity:1" />
    <stop offset="100%" style="stop-color:White;stop-opacity:1" />
  </linearGradient>  
  <g class="branch">
  <polyline points="${treeVisDoc.root.x1} ${treeVisDoc.root.y1}, ${treeVisDoc.root.x2} ${treeVisDoc.root.y2}"/>
  <#list treeVisDoc.branches as branch>
    <polyline <#if branch.highlighted!false>class="highlightedBranch"</#if>
		points="${branch.parentX} ${branch.parentY}, ${branch.cornerX} ${branch.cornerY}, ${branch.childX} ${branch.childY}"/>
  </#list>
  </g>
  <g class="leafNode">
  <#list treeVisDoc.leafNodes as leafNode>
    <#if leafNode.properties.aaValue ??>
     <#if leafNode.properties.multipleResidues??>
	  <g class="multipleAasBackground">
	  <title>Multiple residues: ${leafNode.properties.multipleResidues?split('', 'r')?join(', ')}</title>
	 </#if>
      <circle class="leafNodeAaCircle_${leafNode.properties.aaValue?replace('*', 'STOP')?replace('-', 'GAP')}"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <#assign letterClass = "leafNodeAaText_"+leafNode.properties.aaValue?replace('*', 'STOP')?replace('-', 'GAP')?replace('?', 'X')>
	 <#if leafNode.properties.multipleResidues??>
		  <#assign letterClass = letterClass+" multipleAas">
	 </#if>
     <text class="${letterClass}"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">${leafNode.properties.aaValue}</text>
     <#if leafNode.properties.multipleResidues??>
	  </g>
	 </#if>
	<#elseif leafNode.properties.deletionCode ??>
 	<#if leafNode.properties.deletionCode = "deletion">
      <circle class="leafNodeDeletionCircle_deletion"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <text class="deletionLetter"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">-</text>
	<#elseif leafNode.properties.deletionCode = "no_deletion">
      <circle class="leafNodeDeletionCircle_no_deletion"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <text class="deletionLetter"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">+</text>
	</#if>
	<#elseif leafNode.properties.insertionCode ??>
 	<#if leafNode.properties.insertionCode = "insertion">
      <circle class="leafNodeInsertionCircle_insertion"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <text class="insertionLetter"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">+</text>
	<#elseif leafNode.properties.insertionCode = "insertion_different_nts">
	  <g class="multipleAasBackground">
	  <title>Alternative insertion: ${leafNode.properties.insertedNTs}</title>
      <circle class="leafNodeInsertionCircle_insertion_different_nts"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <text class="insertionLetter"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">?</text>
	  </g>
	<#elseif leafNode.properties.insertionCode = "no_insertion">
      <circle class="leafNodeInsertionCircle_no_insertion"
    	cx="${leafNode.x + leafNodeRadiusPx}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
     <text class="insertionLetter"
	  	x="${leafNode.x + leafNodeRadiusPx}" y="${leafNode.y}" height="${verticalLeafSpacePx}"
	  	dy=".35em">-</text>
	</#if>
	<#else>
      <circle <#if leafNode.highlighted!false>class="leafNodeCircle_highlighted"<#else>class="leafNodeCircle_basic"</#if>
    	cx="${leafNode.x + (leafNodeRadiusPx/2)}" cy="${leafNode.y}" r="${leafNodeRadiusPx/2}"/>
	</#if>

	
  </#list>
  </g>
  <g font-size="${leafTextFontSize}">
  <#list treeVisDoc.leafNodes as leafNode>
    <a id="leafLink${leafNode?index}"><text <#if leafNode.highlighted!false>class="highlightedLeafText"<#else>class="leafText"</#if>
	  	x="${leafNode.x+leafTextGapPx}" y="${leafNode.y}" width="${leafNode.properties.leafTextWidthPx}" height="${verticalLeafSpacePx}"
	  	dy=".35em">${leafNode.properties.leafText?replace('X$', '')}</text></a>
  </#list>
  </g>
  <#if treeVisDoc.collapsedSubtrees??>
  <g class="collapsedSubtree">
  <#list treeVisDoc.collapsedSubtrees as collapsedSubtree>
	<#assign bezierPointXOffsetPx = ((collapsedSubtree.leafX - collapsedSubtree.rootX) * 0.5)>
    <path
		d="M ${collapsedSubtree.rootX} ${collapsedSubtree.rootY}
		C ${collapsedSubtree.rootX + bezierPointXOffsetPx} ${collapsedSubtree.rootY}
		  ${collapsedSubtree.leafX - bezierPointXOffsetPx} ${collapsedSubtree.upperLeafY} 
		  ${collapsedSubtree.leafX} ${collapsedSubtree.upperLeafY} 
		L ${collapsedSubtree.leafX} ${collapsedSubtree.lowerLeafY} 
		C ${collapsedSubtree.leafX - bezierPointXOffsetPx} ${collapsedSubtree.lowerLeafY}
		  ${collapsedSubtree.rootX + bezierPointXOffsetPx} ${collapsedSubtree.rootY}
		  ${collapsedSubtree.rootX} ${collapsedSubtree.rootY} 
		Z"/>
  </#list>
  </g>
  <g font-size="${leafTextFontSize}">
  <#list treeVisDoc.collapsedSubtrees as collapsedSubtree>
	<#if collapsedSubtree.properties.collapsedAlignment ??>
	<a class="collapsedSubtreeLink" href="#project/alignment/${collapsedSubtree.properties.collapsedAlignment}" target="_blank">
	</#if>
    <text class="collapsedSubtreeText"
	  	x="${collapsedSubtree.leafX+leafTextGapPx}" y="${collapsedSubtree.rootY}" width="${collapsedSubtree.properties.collapsedTextWidthPx}" height="${verticalLeafSpacePx}"
	  	dy=".35em">${collapsedSubtree.properties.collapsedText}</text>
	<#if collapsedSubtree.properties.collapsedAlignment ??>
	</a>
	</#if>
  </#list>
  </g>
  </#if>
  <#assign boxTextRCorner = verticalLeafSpacePx / 2>
  <#list treeVisDoc.leafNodes as leafNode>
	<#if leafNode.properties.popupLines ??>
	<#assign boxTextWidth = leafNode.properties.leafPopupTextWidthPx>
	<#assign roundedRectHeight = (leafNode.properties.popupLines?size * verticalLeafSpacePx)+(boxTextRCorner*2)>
	<#assign roundedRectWidth = boxTextWidth+(boxTextRCorner*2)>
	<#assign yCenteringOffset = roundedRectHeight/2>
	<#assign boxSpacer = verticalLeafSpacePx>
	
	<#assign roundedRectX = (leafNode.x+leafTextGapPx+leafNode.properties.leafTextWidthPx+boxSpacer)>
	<#assign roundedRectY = (leafNode.y-yCenteringOffset)>
	<#assign trianglePath = "M"+(leafNode.x+leafTextGapPx+leafNode.properties.leafTextWidthPx+(boxSpacer/2))+" "+leafNode.y+" l"+(1+(boxSpacer/2))+" "+(boxSpacer/2)+" l0 "+(-boxSpacer)+" Z">
	<#if roundedRectX+roundedRectWidth gt svgWidthPx>
		<#assign roundedRectX = (leafNode.x+leafTextGapPx+(leafNode.properties.leafTextWidthPx/2)-(roundedRectWidth/2))>
		<#assign roundedRectY = (leafNode.y-(verticalLeafSpacePx/2 + (boxSpacer*0.75) + roundedRectHeight))>
		<#assign trianglePath = "M"+(leafNode.x+leafTextGapPx+leafNode.properties.leafTextWidthPx/2)+" "+(leafNode.y-(verticalLeafSpacePx/2 + boxSpacer/4))+" l"+(-boxSpacer/2)+" "+(-(1+boxSpacer/2))+" l"+boxSpacer+" 0 Z">
	</#if>
	
	<rect id="popupRect${leafNode?index}" 
		x="${roundedRectX}" 
		y="${roundedRectY}" rx="${boxTextRCorner}" ry="${boxTextRCorner}" 
		fill="#404040" width="${roundedRectWidth}" 
		height="${roundedRectHeight}" 
		visibility="hidden">
	    <set attributeName="visibility" from="hidden" to="visible" begin="leafLink${leafNode?index}.click" end="leafLink${leafNode?index}.mouseout" />
	</rect>
	<path fill="#404040" 
		d="${trianglePath}"
		visibility="hidden">
	    <set attributeName="visibility" from="hidden" to="visible" begin="leafLink${leafNode?index}.click" end="leafLink${leafNode?index}.mouseout"/>
	</path>
	<#list leafNode.properties.popupLines as popupLine>
	<text class="popupText" 
		font-size="${leafTextFontSize}"
		x="${roundedRectX+boxTextRCorner}" 
		y="${roundedRectY+boxTextRCorner+(verticalLeafSpacePx/2)+(popupLine?index * verticalLeafSpacePx)}" 
		width="${boxTextWidth}" dy=".35em" height="${verticalLeafSpacePx}" visibility="hidden">
		${popupLine}
	    <set attributeName="visibility" from="hidden" to="visible" begin="leafLink${leafNode?index}.click" end="leafLink${leafNode?index}.mouseout"/>
	</text>
	</#list>
	</#if>
  </#list>
</svg>
