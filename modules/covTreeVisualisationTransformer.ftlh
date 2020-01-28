<#setting number_format="####0.###">
<#assign treeVisDoc = transformerInput.treeVisualisation>
<#assign svgWidthPx = treeVisDoc.pxWidth>
<#assign svgHeightPx = treeVisDoc.pxHeight>
<#assign verticalLeafSpacePx = treeVisDoc.verticalLeafSpacePx>
<#assign leafTextGapPx = treeVisDoc.leafTextGapPx>
<#assign leafTextFontSize = treeVisDoc.leafTextFontSize>
	
<#assign leafNodeRadiusPx = 2>
	
<svg width="${svgWidthPx}" height="${svgHeightPx}"
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
    <circle <#if leafNode.highlighted!false>class="highlightedLeafNode"</#if>
    	cx="${leafNode.x}" cy="${leafNode.y}" r="${leafNodeRadiusPx}"/>
  </#list>
  </g>
  <g font-size="${leafTextFontSize}">
  <#list treeVisDoc.leafNodes as leafNode>
    <#if leafNode.properties.leafSourceName ?? && leafNode.properties.leafSequenceID ?? >
	<a class="leafLink" href="#project/sequence/${leafNode.properties.leafSourceName}/${leafNode.properties.leafSequenceID}" target="_blank">
	</#if>
    <text <#if leafNode.highlighted!false>class="highlightedLeafText"<#else>class="leafText"</#if>
	  	x="${leafNode.x+leafTextGapPx}" y="${leafNode.y}" width="${leafNode.properties.leafTextWidthPx}" height="${verticalLeafSpacePx}"
	  	dy=".35em">${leafNode.properties.leafText}</text>
    <#if leafNode.properties.leafSourceName ?? && leafNode.properties.leafSequenceID ?? >
	</a>
	</#if>
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



</svg>
