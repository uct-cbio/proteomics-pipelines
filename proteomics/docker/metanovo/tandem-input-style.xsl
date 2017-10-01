<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
                xmlns:GAML="http://www.bioml.com/gaml/" >
<!--
X! tandem default style sheet
Copyright (C) 2003 Ronald C. Beavis
All Rights Reserved
     This source code is distributed under the terms of the
     Artistic License. 
-->

<xsl:template match="/">
  <html>
    <head>
      <title><xsl:value-of select="/bioml/@label" /></title>
      <link rel="stylesheet" href="http://v.thegpm.org/tandem/tandem-style.css" />
    </head>

    <body bgcolor="#FFFFFF">
	<TABLE CELLSPACING="3" CELLPADDING="3">
	<TR>
	<TD WIDTH="500" VALIGN="TOP" ALIGN="LEFT" class="top_note">X! tandem <xsl:value-of select="/bioml/@label" /> method</TD>
	</TR>
	</TABLE>
<table border="0" cellspacing="2" cellpadding="2">
		<xsl:apply-templates select="/bioml/note" />
</table>
    </body>
  </html>
</xsl:template>

<xsl:template match="note">
	<xsl:if test="contains(@type,'input')">
	<xsl:variable name="str_label" select="@label" />
	<xsl:variable name="str_value" select="text()" />
	<TR>
		<TD WIDTH="200" VALIGN="TOP" ALIGN="RIGHT"><I><xsl:value-of select="@label" /></I></TD>
		<TD WIDTH="350" VALIGN="TOP" ALIGN="LEFT"><B><xsl:value-of select="text()" /></B></TD>
	</TR>
	</xsl:if>
	<xsl:if test="contains(@type,'description') ">
	<TR>
		<TD WIDTH="550" VALIGN="TOP" COLSPAN="2" ALIGN="LEFT"><xsl:value-of select="text()" /></TD>
	</TR>
	</xsl:if>
	<xsl:if test="contains(@type,'heading') ">
	<TR>
	<TR>
		<TD WIDTH="200" VALIGN="TOP" ALIGN="CENTER"><B><xsl:value-of select="text()" /></B></TD>
		<TD WIDTH="350" VALIGN="TOP" ALIGN="LEFT"></TD>
	</TR>
	</TR>
	</xsl:if>
</xsl:template>



</xsl:stylesheet>

