/* config - put up track and display configuration page. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "dystring.h"
#include "cheapcgi.h"
#include "htmshell.h"
#include "hdb.h"
#include "hCommon.h"
#include "cart.h"
#include "web.h"
#include "customTrack.h"
#include "hgTracks.h"
#include "hgConfig.h"
#include "jsHelper.h"
#include "imageV2.h"
#include "search.h"
#include "hubConnect.h"
#include "fileUi.h"
#include "trackHub.h"

static void themeDropDown(struct cart* cart)
/* Create drop down for UI themes. 
 * specfied in hg.conf like this
 * browser.theme.modern=background.png,HGStyle
 * */
{
struct slName* themes = cfgNamesWithPrefix("browser.theme.");
//struct slName* themes = cfgNames();
if (themes==NULL)
    return;

hPrintf("<TR><TD>website style:");
hPrintf("<TD style=\"text-align: right\">");

// create labels for drop down box by removing prefix from hg.conf keys
char *labels[50];
struct slName* el;
int i = 0;
el = themes;
for (el = themes; el != NULL && i<50; el = el->next)
    {
    char* name = el->name;
    name = chopPrefix(name); // chop off first two words
    name = chopPrefix(name);
    labels[i] = name;
    i++;
    }

char* currentTheme = cartOptionalString(cart, "theme"); 
hDropList("theme", labels, i, currentTheme);
slFreeList(themes);
hPrintf("</TD>");
}

static void textSizeDropDown()
/* Create drop down for font size. */
{
static char *sizes[] = {"6", "8", "10", "12", "14", "18", "24", "34"};
hDropList(textSizeVar, sizes, ArraySize(sizes), tl.textSize);
}

static void trackConfig(struct track *trackList, struct group *groupList,
	char *groupTarget,  int changeVis)
/* Put up track configurations. If groupTarget is
 * NULL then set visibility for tracks in all groups.  Otherwise,
 * just set it for the given group.  If vis is -2, then visibility is
 * unchanged.  If -1 then set visibility to default, otherwise it should
 * be tvHide, tvDense, etc. */
{
struct group *group;
boolean showedRuler = FALSE;

setRulerMode();
changeTrackVis(groupList, groupTarget, changeVis);

/* Set up ruler mode according to changeVis. */
#ifdef BOB_DOESNT_LIKE
if (changeVis != -2)
    {
    if (groupTarget == NULL ||
        (groupList != NULL && sameString(groupTarget, groupList->name)))
	{
	if (changeVis == -1)
	    rulerMode = tvFull;
	else
	    rulerMode = changeVis;
	}
    }
#endif /* BOB_DOESNT_LIKE */

jsInit();
cgiMakeHiddenVar(configGroupTarget, "none");

// Now all groups are in a single table, divided by an empty borderless row
hPrintf("<TABLE BORDER='0' CELLSPACING='0' class='groupLists'>\n");
for (group = groupList; group != NULL; group = group->next)
    {
    struct trackRef *tr;

    if (group->trackList == NULL)
	continue;

    /* check if group section should be displayed */
    char *otherState;
    char *indicator;
    char *indicatorImg;
    boolean isOpen = !isCollapsedGroup(group);
    collapseGroupGoodies(isOpen, FALSE, &indicatorImg,
                            &indicator, &otherState);
    hPrintf("<TR NOWRAP class='blueToggleBar'>");
    hPrintf("<TH NOWRAP align='left' colspan=3>");
    hPrintf("<table style='width:100%%;'><tr class='noData'><td style='text-align:left;'>");
    hPrintf("\n<A NAME='%sGroup'></A>",group->name);
    hPrintf("<input type=hidden name='%s' id='%s' value=%d>",
            collapseGroupVar(group->name),collapseGroupVar(group->name), (isOpen?0:1));
    hPrintf("<IMG class='toggleButton' "
            "id='%s_togBut' src='%s' alt='%s' title='%s this group'>&nbsp;&nbsp;",
            group->name, indicatorImg, indicator,isOpen?"Collapse":"Expand");
    char jsText[256];
    // TODO XSS filter group->name
    safef(jsText, sizeof jsText, "return vis.toggleForGroup(this,'%s');", group->name);
    char idText[256];
    safef(idText, sizeof idText, "%s_togBut", group->name);
    jsOnEventById("click", idText, jsText);

    hPrintf("<B>&nbsp;%s</B> ", group->label);
    hPrintf("&nbsp;&nbsp;&nbsp;");
    hPrintf("</td><td style='text-align:right;'>\n");
    safef(idText, sizeof idText, "%s_hideAllBut", group->name);
    hPrintf("<INPUT TYPE=SUBMIT NAME=\"%s\" id='%s' VALUE=\"%s\" "
            "title='Hide all tracks in this groups'>",
	    configHideAll, idText, "hide all");
    // TODO XSS filter configGroupTarget
    safef(jsText, sizeof jsText, "document.mainForm.%s.value='%s'; %s",
	    configGroupTarget, group->name, jsSetVerticalPosition("mainForm"));
    jsOnEventById("click", idText, jsText);
    hPrintf(" ");
    safef(idText, sizeof idText, "%s_showAllBut", group->name);
    hPrintf("<INPUT TYPE=SUBMIT NAME=\"%s\" id='%s' VALUE=\"%s\" "
            "title='Show all tracks in this groups'>",
	    configShowAll, idText, "show all");
    jsOnEventById("click", idText, jsText);
    hPrintf(" ");
    safef(idText, sizeof idText, "%s_defaultBut", group->name);
    hPrintf("<INPUT TYPE=SUBMIT NAME=\"%s\" id='%s' VALUE=\"%s\" "
            "title='Show default tracks in this group'>",
	    configDefaultAll, idText, "default");
    jsOnEventById("click", idText, jsText);
    hPrintf(" ");
    /* do not want all the submit buttons named the same.  It is
     * confusing to the javascript submit() function.
     */
    char submitName[256];
    safef(submitName, sizeof(submitName), "%sSubmit", group->name);
    cgiMakeButtonWithMsg(submitName, "submit","Submit your selections and view them in the browser");
    hPrintf("</td></tr></table>\n");
    hPrintf("</TH></TR>\n");

    /* First non-CT, non-hub group gets ruler. */
    if (!showedRuler && !isHubTrack(group->name) &&
                differentString(group->name, "user"))
	{
        showedRuler = TRUE;
        hPrintf("<TR %sid='%s-0'>",(isOpen ? "" : "style='display: none'"), group->name);
	hPrintf("<TD>");
        hPrintf("<A HREF=\"%s?%s=%s&c=%s&g=%s&hgTracksConfigPage=configure\">", hgTrackUiName(),
                cartSessionVarName(), cartSessionId(cart),
                chromName, RULER_TRACK_NAME);
        hPrintf("%s</A>", RULER_TRACK_LABEL);
	hPrintf("</TD><TD>");
	hTvDropDownClass("ruler", rulerMode, FALSE, rulerMode ? "normalText" : "hiddenText");
	hPrintf("</TD><TD>");
	hPrintf("Chromosome position in bases.  (Clicks here zoom in 3x)");
	hPrintf("</TD></TR>\n");
	}
    /* Scan track list to determine which supertracks have visible member
     * tracks, and to insert a track in the list for the supertrack.
     * Sort tracks and supertracks together by priority */
    makeGlobalTrackHash(trackList);
    groupTrackListAddSuper(cart, group);

    if (!withPriorityOverride)
        {
        /* sort hierarchically by priority, considering supertracks */
        struct trackRef *refList = NULL, *ref;
        for (tr = group->trackList; tr != NULL; tr = tr->next)
            {
            struct track *track = tr->track;
            if (tdbIsSuperTrackChild(track->tdb))
                /* ignore supertrack member tracks till supertrack is found */
                continue;
            AllocVar(ref);
            ref->track = track;
            slAddTail(&refList, ref);
            if (tdbIsSuper(track->tdb))
                {
                struct slRef *child = track->tdb->children;
                for (; child != NULL; child=child->next)
                    {
                    struct trackDb *childTdb = child->val;
                    struct track *childTrack = hashFindVal(trackHash, childTdb->track);
                    // Try adding downloadsOnly track
                    if (childTrack == NULL && tdbIsDownloadsOnly(childTdb))
                        {
                        AllocVar(childTrack);           // Fake a track!
                        childTrack->tdb = childTdb;
                        childTrack->hasUi = FALSE;
                        }
                    if (childTrack != NULL)
                        {
                        AllocVar(ref);
                        ref->track = childTrack;
                        slAddTail(&refList, ref);
                        }
                    }
                }
            }
        group->trackList = refList;
        }

    /* Loop through this group and display */
    int rowCount=1;
    for (tr = group->trackList; tr != NULL; tr = tr->next)
        {
        struct track *track = tr->track;
        struct trackDb *tdb = track->tdb;

        hPrintf("<TR %sid='%s-%d'>",(isOpen ? "" : "style='display: none;'"),
                group->name, rowCount++);
        hPrintf("<TD NOWRAP>");
        if (tdbIsSuperTrackChild(tdb))
            /* indent members of a supertrack */
            hPrintf("&nbsp;&nbsp;&nbsp;&nbsp;");

        // Print an icon before the title when one is defined
        hPrintPennantIcon(tdb);

        if (track->hasUi)
            hPrintf("<A TITLE='%s%s...' HREF='%s?%s=%s&g=%s&hgTracksConfigPage=configure'>",
                    tdb->parent ? "Part of super track: " : "Configure ",
                    tdb->parent ? tdb->parent->shortLabel : tdb->shortLabel,
                    hTrackUiForTrack(tdb->track),
                    cartSessionVarName(), cartSessionId(cart), track->track);
        hPrintf(" %s", tdb->shortLabel);
        if (tdbIsSuper(tdb))
            hPrintf("...");
        if (track->hasUi)
	    hPrintf("</A>");
	hPrintf("</TD><TD NOWRAP>");
        if (tdbIsSuperTrackChild(tdb))
            /* indent members of a supertrack */
            hPrintf("&nbsp;&nbsp;&nbsp;&nbsp;");

	/* If track is not on this chrom print an informational
	   message for the user. */
        if (tdbIsDownloadsOnly(tdb))    // No vis display for downloadsOnly
            hPrintf("<A TITLE='Downloadable files...' HREF='%s?%s=%s&g=%s'>Downloads</A>",
                    hgFileUiName(),cartSessionVarName(), cartSessionId(cart), tdb->track);
        else if (hTrackOnChrom(track->tdb, chromName))
	    {
            if (tdbIsSuper(track->tdb))
                {
                /* supertrack dropdown is hide/show */
                superTrackDropDown(cart, track->tdb, 1);
                }
            else
                {
                /* check for option of limiting visibility to one mode */
                hTvDropDownClassVisOnly(track->track, track->visibility,
                                        rTdbTreeCanPack(track->tdb),
                                        (track->visibility == tvHide) ? "hiddenText" : "normalText",
                                        trackDbSetting(track->tdb, "onlyVisibility"));
                }
	    }
        else
	    hPrintf("[No data-%s]", chromName);
	hPrintf("</TD><TD NOWRAP>");
        hPrintf("%s", tdb->longLabel);
	hPrintf("</TD></TR>\n");
	}
    hPrintf("<tr class='noData'><td colspan=3>");
    cgiDown(0.9);
    hPrintf("</td></tr>\n");
    }
hPrintf("</TABLE>\n");
}

static int addDownloadOnlyTracks(char *db,struct group **pGroupList,struct track **pTrackList)
// Download only tracks are not normaly incorporated into the grou and track lists
{
int count = 0;
struct track *track = NULL;
struct group *group = NULL;
struct trackDb *tdbList = hTrackDb(db);
struct trackDb *tdb = tdbList;
for (;tdb != NULL; tdb = tdb->next)
    {
    if (!tdbIsDownloadsOnly(tdb)
    || tdbIsFolderContent(tdb)
    || tdbIsCompositeChild(tdb))
        continue;

    // Must find group
    if (tdb->grp == NULL)
        continue;

    for (group = *pGroupList;group != NULL; group = group->next)
        {
        if (sameWord(group->name,tdb->grp))
            break;
        }
    if (group == NULL)
        continue;

    // Make the track
    track = trackFromTrackDb(tdb);
    track->group = group;
    track->groupName = cloneString(group->name);
    slAddHead(pTrackList,track);
    count++;
    }

if (count > 0)
    {
    // Going to have to make all new group->trackLists
    slSort(pGroupList, gCmpPriority);
    for (group = *pGroupList;group != NULL; group = group->next)
        slFreeList(&group->trackList);

    // Sort the tracks anew and add each on into it's group.
    slSort(pTrackList, tgCmpPriority);
    for (track = *pTrackList; track != NULL; track = track->next)
        {
        struct trackRef *tr;
        AllocVar(tr);
        tr->track = track;
        slAddHead(&track->group->trackList, tr);
        }

    /* Straighten things out, clean up, and go home. */
    for (group = *pGroupList;group != NULL; group = group->next)
        slReverse(&group->trackList);
    }
return count;
}


void configInitTrackList(
    int vis, 
    char **pGroupTarget,
    struct track **pTrackList,
    struct track **pIdeoTrack,
    struct group **pGroupList
)
{
char *groupTarget = NULL;
struct track *trackList =  NULL;
struct track *ideoTrack = NULL;
struct group *groupList = NULL;

withPriorityOverride = cartUsualBoolean(cart, configPriorityOverride, FALSE);

/* Get track list and group them. */
ctList = customTracksParseCart(database, cart, &browserLines, &ctFileName);
trackList = getTrackList(&groupList, vis);

if (trackHash == NULL)
    trackHash = makeGlobalTrackHash(trackList);
// Subtrack settings must be removed when composite/view settings are updated
parentChildCartCleanup(trackList,cart,oldVars);

addDownloadOnlyTracks(database,&groupList,&trackList);

/* The ideogram for some reason is considered a track.
 * We don't really want to process it as one though, so
 * we see if it's there, and if necessary remove it. */
ideoTrack = chromIdeoTrack(trackList);
if (ideoTrack != NULL)
    removeTrackFromGroup(ideoTrack);

/* Fetch group to change on if any from CGI,
 * and remove var so it doesn't get used again. */
groupTarget = cloneString(cartUsualString(cart, configGroupTarget, ""));
cartRemove(cart, configGroupTarget);
if (sameString(groupTarget, "none"))
    freez(&groupTarget);

*pGroupTarget = groupTarget;
*pTrackList = trackList;
*pIdeoTrack = ideoTrack;
*pGroupList = groupList;

}

void configPageSetTrackVis(int vis)
/* Do config page after setting track visibility. If vis is -2, then visibility
 * is unchanged.  If -1 then set visibility to default, otherwise it should
 * be tvHide, tvDense, etc. */
{
char *groupTarget;
struct track *trackList;
struct track *ideoTrack;
struct group *groupList;

configInitTrackList(vis, &groupTarget, &trackList, &ideoTrack, &groupList);

struct dyString *title = dyStringNew(0);

dyStringPrintf(title, "Configure Image");

hPrintf("<FORM ACTION=\"%s\" NAME=\"mainForm\" METHOD=%s>\n", hgTracksName(),
	cartUsualString(cart, "formMethod", "POST"));
webStartWrapperDetailedNoArgs(cart, database, "", title->string, FALSE, FALSE, FALSE, FALSE);
cartSaveSession(cart);

hPrintf("<INPUT TYPE=HIDDEN NAME=\"hgTracksConfigPage\" VALUE=\"\">");
/* do not want all the submit buttons named the same thing, this one is: */
cgiMakeButton("topSubmit", "submit");

// 3 column table
hPrintf("<TABLE style=\"border:0px; \">\n");
hPrintf("<TR><TD>image width:");
hPrintf("<TD style=\"text-align: right\">");
hIntVar("pix", tl.picWidth, 4);
hPrintf("<TD>pixels</TR>");

hPrintf("<TR><TD>label area width:");
hPrintf("<TD style=\"text-align: right\">");
hIntVar("hgt.labelWidth", leftLabelWidthChars, 2);
hPrintf("<TD>characters<TD></TR>");

hPrintf("<TR><TD>text size:");
hPrintf("<TD style=\"text-align: right\">");
textSizeDropDown();
hPrintf("</TD>");
if (trackLayoutInclFontExtras())
    {
    hPrintf("<TD>");
    char *defaultStyle = cartUsualString(cart, "fontType", "medium");
    cartMakeRadioButton(cart, "fontType", "medium", defaultStyle);
    hPrintf("&nbsp;medium&nbsp;");
    cartMakeRadioButton(cart, "fontType", "fixed", defaultStyle);
    hPrintf("&nbsp;fixed&nbsp;");
    cartMakeRadioButton(cart, "fontType", "bold", defaultStyle);
    hPrintf("&nbsp;bold&nbsp;");
    hPrintf("&nbsp;");
    hPrintf("</TD>");
    }
hPrintf("</TR>");

themeDropDown(cart);

hTableStart();
if (ideoTrack != NULL)
    {
    hPrintf("<TR><TD>");
    hCheckBox("ideogram", cartUsualBoolean(cart, "ideogram", TRUE));
    hPrintf("</TD><TD>");
    hPrintf("Display chromosome ideogram above main graphic");
    hPrintf("</TD></TR>\n");
    }
hPrintf("<TR><TD>");
hCheckBox("guidelines", cartUsualBoolean(cart, "guidelines", TRUE));
hPrintf("</TD><TD>");
hPrintf("Show light blue vertical guidelines, or light red vertical  window separators in multi-region view");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("leftLabels", cartUsualBoolean(cart, "leftLabels", TRUE));
hPrintf("</TD><TD>");
hPrintf("Display labels to the left of items in tracks");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("centerLabels", cartUsualBoolean(cart, "centerLabels", TRUE));
hPrintf("</TD><TD>");
hPrintf("Display description above each track");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("trackControlsOnMain", cartUsualBoolean(cart, "trackControlsOnMain", TRUE));
hPrintf("</TD><TD>");
hPrintf("Show track controls under main graphic");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("nextItemArrows", cartUsualBoolean(cart, "nextItemArrows", FALSE));
hPrintf("</TD><TD>");
hPrintf("Next/previous item navigation");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("nextExonArrows", cartUsualBoolean(cart, "nextExonArrows", TRUE));
hPrintf("</TD><TD>");
hPrintf("Next/previous exon navigation");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("exonNumbers", cartUsualBoolean(cart, "exonNumbers", TRUE));
hPrintf("</TD><TD>");
hPrintf("Show exon numbers");
hPrintf("</TD></TR>\n");

hPrintf("<TR><TD>");
hCheckBox("enableHighlightingDialog", cartUsualBoolean(cart, "enableHighlightingDialog", TRUE));
hPrintf("</TD><TD>");
hPrintf("Enable highlight with drag-and-select "
        "(if unchecked, drag-and-select always zooms to selection)");
hPrintf("</TD></TR>\n");
hTableEnd();


cgiDown(0.9);

char buf[256];
char *freeze = hFreezeFromDb(database);
if (freeze == NULL)
    safef(buf, sizeof buf, "Configure Tracks on %s %s: %s",
	  organization, browserName, trackHubSkipHubName(organism));
else if (stringIn(database, freeze))
    safef(buf, sizeof buf, "Configure Tracks on %s %s: %s %s",
	  organization, browserName, trackHubSkipHubName(organism), freeze);
else
    safef(buf, sizeof buf, "Configure Tracks on %s %s: %s %s (%s)",
	  organization, browserName, trackHubSkipHubName(organism), freeze, trackHubSkipHubName(database));
webNewSection("%s",buf);
hPrintf("Tracks: ");
if (isSearchTracksSupported(database,cart))
    {
    cgiMakeButtonWithMsg(TRACK_SEARCH, TRACK_SEARCH_BUTTON,TRACK_SEARCH_HINT);
    hPrintf(" ");
    }
cgiMakeButtonWithMsg(configHideAll, "hide all","Hide all tracks in this genome assembly");
hPrintf(" ");
cgiMakeButtonWithMsg(configShowAll, "show all","Show all tracks in this genome assembly");
hPrintf(" ");
cgiMakeButtonWithMsg(configDefaultAll, "default","Display only default tracks");
hPrintf("&nbsp;&nbsp;&nbsp;Groups:  ");
hButtonWithOnClick("hgt.collapseGroups", "collapse all", "Collapse all track groups",
                   "return vis.expandAllGroups(false)");
hPrintf(" ");
hButtonWithOnClick("hgt.expandGroups", "expand all", "Expand all track groups",
                   "return vis.expandAllGroups(true)");
hPrintf("<div style='margin-top:.2em; margin-bottom:.9em;'>Control track and group visibility "
        "more selectively below.</div>");
trackConfig(trackList, groupList, groupTarget, vis);

dyStringFree(&title);
freez(&groupTarget);
webEndSectionTables();
hPrintf("</FORM>");
}

void configPage()
/* Put up configuration page. */
{
configPageSetTrackVis(-2);
}

// TODO GALT there is duplication still between config and configMultiRegionPageSetTrackVis
//  that could maybe be addressed by pulling the code that initializes the tracklist,
//  and the code that draws the multi-region options, into 2 functions to be called by each.

void configMultiRegionPage()
/* Do multi-region config page after setting track visibility. If vis is -2, then visibility
 * is unchanged.  If -1 then set visibility to default, otherwise it should
 * be tvHide, tvDense, etc. */
{
char *groupTarget;
struct track *trackList;
struct track *ideoTrack;
struct group *groupList;
int vis = -2;

configInitTrackList(vis, &groupTarget, &trackList, &ideoTrack, &groupList);

hPrintf("<FORM ACTION=\"%s\" NAME=\"mainForm\" METHOD=%s>\n", hgTracksName(),
	cartUsualString(cart, "formMethod", "POST"));

webStartWrapperDetailedNoArgs(cart, database, "", "", FALSE, FALSE, FALSE, FALSE);

cartSaveSession(cart);


hPrintf("<BR>\n");

hTableStart();

virtModeType = cartUsualString(cart, "virtModeType", virtModeType);

hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "default", sameWord("default", virtModeType));
hPrintf("</TD><TD>");
hPrintf("Exit multi-region mode");
hPrintf("</TD></TR>\n");

struct sqlConnection *conn = NULL;
if (!trackHubDatabase(database))  // no db conn for assembly hubs 
    conn = hAllocConn(database);

// Do we have a gene table for exonMostly?
findBestEMGeneTable(trackList);
if (emGeneTable)
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "exonMostly", sameWord("exonMostly", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("Show exons using %s. &nbsp;&nbsp; Use padding of: ", emGeneTrack->shortLabel);
    hIntVar("emPadding", cartUsualInt(cart, "emPadding", emPadding), 3);
    hPrintf(" bases.");
    hPrintf("</TD></TR>\n");
    }

if (emGeneTable)
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "geneMostly", sameWord("geneMostly", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("Show genes using %s. &nbsp;&nbsp; Use padding of: ", emGeneTrack->shortLabel);
    hIntVar("gmPadding", cartUsualInt(cart, "gmPadding", gmPadding), 3);
    hPrintf(" bases.");
    hPrintf("</TD></TR>\n");
    }

/* obsolete    
if (conn && sqlTableExists(conn,"knownCanonical"))
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "kcGenes", sameWord("kcGenes", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("Show gene regions genome-wide.");
    hPrintf("</TD></TR>\n");
    }
*/

hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "customUrl", sameWord("customUrl", virtModeType));
hPrintf("</TD><TD>");
hPrintf("Custom regions from BED URL: ");
hTextVar("multiRegionsBedUrl", cartUsualString(cart, "multiRegionsBedUrl", multiRegionsBedUrl), 60);
hPrintf("</TD></TR>\n");


/* The AllChroms option will be released in future
if (emGeneTable && sqlTableExists(conn, emGeneTable))
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "singleTrans", sameWord("singleTrans", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("Show only one transcript using an ID from %s : ", emGeneTrack->shortLabel);
    char *trans = cartUsualString(cart, "singleTransId", singleTransId);
    char sql[1024];
    sqlSafef(sql, sizeof sql, "select name from %s where name='%s'", emGeneTable, trans);
    char *result = sqlQuickString(conn, sql);
    if (!result)
	{
	sqlSafef(sql, sizeof sql, "select name from %s limit 1", emGeneTable);
	trans = sqlQuickString(conn, sql);
	}
    hTextVar("singleTransId", trans, 20);
    hPrintf("</TD></TR>\n");
    }
*/

if (conn && sqlTableExists(conn, "altLocations"))
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "singleAltHaplo", sameWord("singleAltHaplo", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("Show one alternate haplotype, placed on its chromosome, using ID: ");
    char *haplo = cartUsualString(cart, "singleAltHaploId", singleAltHaploId);
    char sql[1024];
    sqlSafef(sql, sizeof sql, "select name from altLocations where name='%s'", haplo);
    char *result = sqlQuickString(conn, sql);
    if (!result)
	{
	sqlSafef(sql, sizeof sql, "select name from altLocations limit 1");
	haplo = sqlQuickString(conn, sql);
	}
    hTextVar("singleAltHaploId", haplo, 20);
    hPrintf("</TD></TR>\n");
    }

/* disable demo for now
if (sameString(database,"hg19") || sameString(database, "hg38"))
    {
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "demo1", sameWord("demo1", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("demo1 two windows on two chroms (default pos on chr21, and same loc on chr22)");
    hPrintf("</TD></TR>\n");
    }
*/


/* Disabled for now 
hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "demo2", sameWord("demo2", virtModeType));
hPrintf("</TD><TD>");
hPrintf("demo2 multiple "); 
hIntVar("demo2NumWindows", cartUsualInt(cart, "demo2NumWindows", demo2NumWindows), 3);
hPrintf(" windows on one chrom chr21 def posn, window size ");
hIntVar("demo2WindowSize", cartUsualInt(cart, "demo2WindowSize", demo2WindowSize), 3);
hPrintf(" and step size ");
hIntVar("demo2StepSize", cartUsualInt(cart, "demo2StepSize", demo2StepSize), 3);
hPrintf(" exon-like");
hPrintf("</TD></TR>\n");
*/

/* The AllChroms option will be released in future
if (conn)  // requires chromInfo from database. 
    { // TODO allow it to use assembly hubs via trackHubAllChromInfo() ?
    hPrintf("<TR><TD>");
    cgiMakeRadioButton("virtModeType", "allChroms", sameWord("allChroms", virtModeType));
    hPrintf("</TD><TD>");
    hPrintf("<br>Show all chromosomes.<br><span style='color:red'>Warning:</span> Turn off all tracks except bigBed, bigWig, and very sparse tracks.<br>Press Hide All to hide all tracks.");
    hPrintf("</TD></TR>\n");
    }
*/


/* Disabled for now 
hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "demo4", sameWord("demo4", virtModeType));
hPrintf("</TD><TD>");
hPrintf("demo4 multiple (311) windows showing exons from TITIN gene uc031rqd.1.");
hPrintf("</TD></TR>\n");
*/

/* Disabled for now 
hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "demo5", sameWord("demo5", virtModeType));
hPrintf("</TD><TD>");
hPrintf("demo5 alt locus on hg38. Shows alt chrom surrounded by regions of same size from reference genome.");
hPrintf("</TD></TR>\n");
*/

/* Disabled for now 
hPrintf("<TR><TD>");
cgiMakeRadioButton("virtModeType", "demo6", sameWord("demo6", virtModeType));
hPrintf("</TD><TD>");
hPrintf("demo6 shows zoomed in exon-exon junction from SOD1 gene, between exon1 and exon2.");
hPrintf("</TD></TR>\n");
*/


hTableEnd();

hPrintf("<BR>\n");
hPrintf("<TABLE style=\"border:0px; \">\n");
hPrintf("<TR><TD>");
hCheckBox("emAltHighlight", cartUsualBoolean(cart, "emAltHighlight", FALSE));
hPrintf("</TD><TD>");
hPrintf("Highlight alternating regions in multi-region view");
hPrintf("</TD></TR>\n");
hPrintf("</TABLE>\n");

hPrintf("<BR>\n");
hPrintf("<TABLE style=\"border:0px;width:650px \">\n");
hPrintf("<TR><TD>");
cgiMakeButton("topSubmit", "submit");
hPrintf("</TD><TD align=right>");
hPrintf("<A HREF=\"../goldenPath/help/multiRegionHelp.html\" target=_blank>Help</A>\n");
hPrintf("</TD></TR>\n");
hPrintf("</TABLE>\n");

hFreeConn(&conn);

cgiDown(0.9);

freez(&groupTarget);
webEndSectionTables();
hPrintf("</FORM>");
}

