<html>
<head>
<meta charset="UTF-8">
<style>
${getResourceAsString('covPrimerProbeReportTransformer.css')?no_esc}
</style>
</head>
<body>

<div id="logo">
	<img id="cvrLogo" alt="CVR logo" src="data:image/jpeg;base64,${getResourceAsBase64('cvrLogo.png')}" />
</div>

<div id="title">
<h2>SARS-CoV-2 Primer/probe mismatch report</h2>
</div>

<div id="overview">
<table>
<colgroup>
    <col class="tableHeader"/>
    <col/>
</colgroup>

<tr><td><b>File path</b></td><td>${covPPReport.fastaFilePath}</td></tr>
<tr><td><b>Query sequence ID</b></td><td>${covPPReport.sequenceID}</td></tr>
<tr><td><b>Report generation date</b></td><td>${covPPReport.reportGenerationDate}</td></tr>
<tr><td><b>GLUE engine version</b></td><td>${covPPReport.glueVersion}</td></tr>
<tr><td><b>CoV-GLUE project version</b></td><td>${covPPReport.projectVersion}</td></tr>
<tr><td><b>Reporting software author</b></td><td>Josh Singer &lt;josh.singer@glasgow.ac.uk&gt;</td></tr>
</table>
</div>

<p>
<h4>Assay reports</h4>
<table id="assayReports">
<thead>
<col width="20%"></col>
<col width="15%"></col>
<col width="15%"></col>
<col width="20%"></col>
<col width="10%"></col>
<col width="20%"></col>

<tr class="tableHeader">
  <td><b>Publication</b></td>
  <td><b>Assay</b></td>
  <td><b>Primer/probe</b></td>
  <td><b>Primer/probe sequence</b></td>
  <td><b>Location on reference</b></td>
  <td><b>Query sequence issues</b></td>
</tr>
</thead>
<#assign knownIssues = false>
<tbody>
<#list covPPReport.assays as assay>
<tr>
<#if assay.primersWithIssues gt 0>
<td rowspan="${assay.primersWithIssues}">
<#else>
<td>
</#if>
<a href="${assay.url}" target="_blank">${assay.organisation}</a>
</td>
<#if assay.primersWithIssues gt 0>
<td rowspan="${assay.primersWithIssues}">
<#else>
<td>
</#if>
${assay.display_name}</td>
<#if assay.primersWithIssues = 0>
<td colspan="4">No issues detected</td>
</tr>
<#else>
<#list assay.ppObjs as ppObj>
<#assign firstIssueFound = false>
<#if ppObj.numKnownIssues gt 0>
<#assign knownIssues = true>
</#if>
<#if ppObj.numIssues gt 0>
<#if firstIssueFound>
<tr>
<#else>
<#assign firstIssueFound = true>
</#if>
<#assign class = "issue">
<#if ppObj.numIssues = ppObj.numKnownIssues>
<#assign class = "">
</#if>
<td class="${class}">${ppObj.display_name}</td>
<td class="${class}">${ppObj.sequence_to_scan}</td>
<td class="${class}">${ppObj.ref_start?c}-${ppObj.ref_end?c}</td>
<td class="${class}">
<#list ppObj.displayIssues as displayIssue>
<#if displayIssue?index gt 0>
<br/>
</#if>
<span>${displayIssue}</span>
</#list>
</td>
</tr>
</#if>
</#list>
</#if>
</#list>
</tbody>
</table>
<#if knownIssues>
<p>
<small>* This is a known issue relating to the primer/probe design rather than a sequence polymorphism.</small>
</p>
</#if>
</body>
</html>