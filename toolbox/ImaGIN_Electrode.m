<!DOCTYPE html>
<html lang='en'>
<head prefix='og: http://ogp.me/ns#'>
<meta charset='utf-8'>
<meta content='IE=edge' http-equiv='X-UA-Compatible'>
<meta content='object' property='og:type'>
<meta content='GitLab' property='og:site_name'>
<meta content='toolbox/ImaGIN_Electrode.m · c4081fa6f23394f3f85ab37ca70e4f12c0752677 · ftract_dev / ImaGIN2' property='og:title'>
<meta content='Reorganization of ImaGIN code (future SPM toolbox)' property='og:description'>
<meta content='https://gin11-git.ujf-grenoble.fr/assets/gitlab_logo-cdf021b35c4e6bb149e26460f26fae81e80552bc879179dd80e9e9266b14e894.png' property='og:image'>
<meta content='https://gin11-git.ujf-grenoble.fr/ftract_dev/ImaGIN2/blob/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m' property='og:url'>
<meta content='summary' property='twitter:card'>
<meta content='toolbox/ImaGIN_Electrode.m · c4081fa6f23394f3f85ab37ca70e4f12c0752677 · ftract_dev / ImaGIN2' property='twitter:title'>
<meta content='Reorganization of ImaGIN code (future SPM toolbox)' property='twitter:description'>
<meta content='https://gin11-git.ujf-grenoble.fr/assets/gitlab_logo-cdf021b35c4e6bb149e26460f26fae81e80552bc879179dd80e9e9266b14e894.png' property='twitter:image'>

<title>toolbox/ImaGIN_Electrode.m · c4081fa6f23394f3f85ab37ca70e4f12c0752677 · ftract_dev / ImaGIN2 · GitLab</title>
<meta content='Reorganization of ImaGIN code (future SPM toolbox)' name='description'>
<link rel="shortcut icon" type="image/x-icon" href="/assets/favicon-075eba76312e8421991a0c1f89a89ee81678bcde72319dd3e8047e2a47cd3a42.ico" />
<link rel="stylesheet" media="all" href="/assets/application-6867d81d1da4ea60f24c4e4f0335d21126326e94418c9de357cf3c397926ee03.css" />
<link rel="stylesheet" media="print" href="/assets/print-68eed6d8135d858318821e790e25da27b2b4b9b8dbb1993fa6765d8e2e3e16ee.css" />
<script src="/assets/application-95adfce412d4f6bc8ed3f7894a4e1572f04202ad6090f05b74ccf98757f61889.js"></script>
<meta name="csrf-param" content="authenticity_token" />
<meta name="csrf-token" content="MnV3pnc/1tyDyyGEdLqz85rlVeonM4UKQebpshNH3YQ6sr65jzOJ/gPS/pMIQ4kWYhSp5VGMt0y0W0s9Cr1y8A==" />
<script>
//<![CDATA[
window.gon={};gon.api_version="v3";gon.default_avatar_url="https://gin11-git.ujf-grenoble.fr/assets/no_avatar-07eeb128b993e74003e8664cff0b8e1e7234cec0443766a6763df32ca3472c23.png";gon.default_issues_tracker="gitlab";gon.max_file_size=10;gon.relative_url_root="";gon.shortcuts_path="/help/shortcuts";gon.user_color_scheme="white";gon.current_user_id=48;gon.api_token="De1V__LWH_G7kZcreDyp";
//]]>
</script>
<meta content='origin-when-cross-origin' name='referrer'>
<meta content='width=device-width, initial-scale=1, maximum-scale=1' name='viewport'>
<meta content='#474D57' name='theme-color'>
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-2d64ecc33893baab71adc15ff19a803a59911cc2651fb9d4d62af1379ff89aab.png" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-d08897d57e1bc9400024ac15465340e832a8e7b166b91624138d48ea2c739596.png" sizes="76x76" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-retina-81446c57f3351d1dacd0fb5f23ced74ba63d3878810bedea343999c6a12b3915.png" sizes="120x120" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-retina-e2a776da039936ac240e76615add47b25ab77c75a5fa2d1b3907f83d5502b911.png" sizes="152x152" />
<link color='rgb(226, 67, 41)' href='/assets/logo-d36b5212042cebc89b96df4bf6ac24e43db316143e89926c0db839ff694d2de4.svg' rel='mask-icon'>
<meta content='/assets/msapplication-tile-49c9c46afd2ab9bbf35e69138bc62f8c31fa55584bd494956ac6e58e6aadc813.png' name='msapplication-TileImage'>
<meta content='#30353E' name='msapplication-TileColor'>




<style>
  [data-user-is] {
    display: none !important;
  }
  
  [data-user-is="48"] {
    display: block !important;
  }
  
  [data-user-is="48"][data-display="inline"] {
    display: inline !important;
  }
  
  [data-user-is-not] {
    display: block !important;
  }
  
  [data-user-is-not][data-display="inline"] {
    display: inline !important;
  }
  
  [data-user-is-not="48"] {
    display: none !important;
  }
</style>

</head>

<body class='ui_charcoal' data-page='projects:blob:show'>
<script>
  window.project_uploads_path = "/ftract_dev/ImaGIN2/uploads";
  window.markdown_preview_path = "/ftract_dev/ImaGIN2/markdown_preview";
</script>

<header class='header-expanded navbar navbar-fixed-top navbar-gitlab'>
<div class='container-fluid'>
<div class='header-content'>
<button class='side-nav-toggle' type='button'>
<span class='sr-only'>Toggle navigation</span>
<i class="fa fa-bars"></i>
</button>
<button class='navbar-toggle' type='button'>
<span class='sr-only'>Toggle navigation</span>
<i class="fa fa-angle-left"></i>
</button>
<div class='navbar-collapse collapse'>
<ul class='nav navbar-nav'>
<li class='hidden-sm hidden-xs'>
<div class='has-location-badge search search-form'>
<form class="navbar-form" action="/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
<div class='search-input-container'>
<div class='search-location-badge'>
<span class='location-badge'>
<i class='location-text'>
This project
</i>
</span>
</div>
<div class='search-input-wrap'>
<div class='dropdown' data-url='/search/autocomplete'>
<input type="search" name="search" id="search" placeholder="Search" class="search-input dropdown-menu-toggle" spellcheck="false" tabindex="1" autocomplete="off" data-toggle="dropdown" />
<div class='dropdown-menu dropdown-select'>
<div class="dropdown-content"><ul>
<li>
<a class='is-focused dropdown-menu-empty-link'>
Loading...
</a>
</li>
</ul>
</div><div class="dropdown-loading"><i class="fa fa-spinner fa-spin"></i></div>
</div>
<i class='search-icon'></i>
<i class='clear-icon js-clear-input'></i>
</div>
</div>
</div>
<input type="hidden" name="group_id" id="group_id" />
<input type="hidden" name="project_id" id="search_project_id" value="30" />
<input type="hidden" name="search_code" id="search_code" value="true" />
<input type="hidden" name="repository_ref" id="repository_ref" value="c4081fa6f23394f3f85ab37ca70e4f12c0752677" />

<div class='search-autocomplete-opts hide' data-autocomplete-path='/search/autocomplete' data-autocomplete-project-id='30' data-autocomplete-project-ref='c4081fa6f23394f3f85ab37ca70e4f12c0752677'></div>
</form>

</div>

</li>
<li class='visible-sm visible-xs'>
<a title="Search" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/search"><i class="fa fa-search"></i>
</a></li>
<li>
<a title="Todos" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/dashboard/todos"><i class="fa fa-bell fa-fw"></i>
<span class='badge todos-pending-count'>
13
</span>
</a></li>
<li>
<a title="New project" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/projects/new"><i class="fa fa-plus fa-fw"></i>
</a></li>
<li>
<a class="logout" title="Sign out" data-toggle="tooltip" data-placement="bottom" data-container="body" rel="nofollow" data-method="delete" href="/users/sign_out"><i class="fa fa-sign-out"></i>
</a></li>
</ul>
</div>
<h1 class='title'><a href="/groups/ftract_dev">ftract_dev</a> / <a class="project-item-select-holder" href="/ftract_dev/ImaGIN2">ImaGIN2</a><i data-target=".js-dropdown-menu-projects" data-toggle="dropdown" class="fa fa-chevron-down dropdown-toggle-caret js-projects-dropdown-toggle"></i> &middot; <a href="/ftract_dev/ImaGIN2/tree/c4081fa6f23394f3f85ab37ca70e4f12c0752677">Files</a></h1>
<div class='js-dropdown-menu-projects'>
<div class='dropdown-menu dropdown-select dropdown-menu-projects'>
<div class="dropdown-title"><span>Go to a project</span><button class="dropdown-title-button dropdown-menu-close" aria-label="Close" type="button"><i class="fa fa-times dropdown-menu-close-icon"></i></button></div>
<div class="dropdown-input"><input type="search" id="" class="dropdown-input-field" placeholder="Search your projects" /><i class="fa fa-search dropdown-input-search"></i><i role="button" class="fa fa-times dropdown-input-clear js-dropdown-input-clear"></i></div>
<div class="dropdown-content"></div>
<div class="dropdown-loading"><i class="fa fa-spinner fa-spin"></i></div>
</div>
</div>

</div>
</div>
</header>

<script>
  var findFileURL = "/ftract_dev/ImaGIN2/find_file/c4081fa6f23394f3f85ab37ca70e4f12c0752677";
</script>

<div class='page-sidebar-expanded page-with-sidebar'>
<div class='nicescroll sidebar-expanded sidebar-wrapper'>
<div class='header-logo'>
<a id='logo'>
<svg width="36" height="36" id="tanuki-logo">
  <path id="tanuki-right-ear" class="tanuki-shape" fill="#e24329" d="M2 14l9.38 9v-9l-4-12.28c-.205-.632-1.176-.632-1.38 0z"/>
  <path id="tanuki-left-ear" class="tanuki-shape" fill="#e24329" d="M34 14l-9.38 9v-9l4-12.28c.205-.632 1.176-.632 1.38 0z"/>
  <path id="tanuki-nose" class="tanuki-shape" fill="#e24329" d="M18,34.38 3,14 33,14 Z"/>
  <path id="tanuki-right-eye" class="tanuki-shape" fill="#fc6d26" d="M18,34.38 11.38,14 2,14 6,25Z"/>
  <path id="tanuki-left-eye" class="tanuki-shape" fill="#fc6d26" d="M18,34.38 24.62,14 34,14 30,25Z"/>
  <path id="tanuki-right-cheek" class="tanuki-shape" fill="#fca326" d="M2 14L.1 20.16c-.18.565 0 1.2.5 1.56l17.42 12.66z"/>
  <path id="tanuki-left-cheek" class="tanuki-shape" fill="#fca326" d="M34 14l1.9 6.16c.18.565 0 1.2-.5 1.56L18 34.38z"/>
</svg>

</a>
<a class="gitlab-text-container-link" title="Dashboard" id="js-shortcuts-home" href="/"><div class='gitlab-text-container'>
<h3>GitLab</h3>
</div>
</a></div>
<ul class='nav nav-sidebar'>
<li class=""><a title="Go to group" class="back-link" href="/groups/ftract_dev"><i class="fa fa-caret-square-o-left fa-fw"></i>
<span>
Go to group
</span>
</a></li><li class='separate-item'></li>
<li class="home"><a title="Project" class="shortcuts-project" href="/ftract_dev/ImaGIN2"><i class="fa fa-bookmark fa-fw"></i>
<span>
Project
</span>
</a></li><li class=""><a title="Activity" class="shortcuts-project-activity" href="/ftract_dev/ImaGIN2/activity"><i class="fa fa-dashboard fa-fw"></i>
<span>
Activity
</span>
</a></li><li class="active"><a title="Files" class="shortcuts-tree" href="/ftract_dev/ImaGIN2/tree/c4081fa6f23394f3f85ab37ca70e4f12c0752677"><i class="fa fa-files-o fa-fw"></i>
<span>
Files
</span>
</a></li><li class=""><a title="Commits" class="shortcuts-commits" href="/ftract_dev/ImaGIN2/commits/c4081fa6f23394f3f85ab37ca70e4f12c0752677"><i class="fa fa-history fa-fw"></i>
<span>
Commits
</span>
</a></li><li class=""><a title="Pipelines" class="shortcuts-pipelines" href="/ftract_dev/ImaGIN2/pipelines"><i class="fa fa-ship fa-fw"></i>
<span>
Pipelines
<span class='count ci_counter'>0</span>
</span>
</a></li><li class=""><a title="Builds" class="shortcuts-builds" href="/ftract_dev/ImaGIN2/builds"><i class="fa fa-cubes fa-fw"></i>
<span>
Builds
<span class='count builds_counter'>0</span>
</span>
</a></li><li class=""><a title="Graphs" class="shortcuts-graphs" href="/ftract_dev/ImaGIN2/graphs/c4081fa6f23394f3f85ab37ca70e4f12c0752677"><i class="fa fa-area-chart fa-fw"></i>
<span>
Graphs
</span>
</a></li><li class=""><a title="Milestones" href="/ftract_dev/ImaGIN2/milestones"><i class="fa fa-clock-o fa-fw"></i>
<span>
Milestones
</span>
</a></li><li class=""><a title="Issues" class="shortcuts-issues" href="/ftract_dev/ImaGIN2/issues"><i class="fa fa-exclamation-circle fa-fw"></i>
<span>
Issues
<span class='count issue_counter'>7</span>
</span>
</a></li><li class=""><a title="Merge Requests" class="shortcuts-merge_requests" href="/ftract_dev/ImaGIN2/merge_requests"><i class="fa fa-tasks fa-fw"></i>
<span>
Merge Requests
<span class='count merge_counter'>0</span>
</span>
</a></li><li class=""><a title="Members" class="team-tab tab" href="/ftract_dev/ImaGIN2/project_members"><i class="fa fa-users fa-fw"></i>
<span>
Members
</span>
</a></li><li class=""><a title="Labels" href="/ftract_dev/ImaGIN2/labels"><i class="fa fa-tags fa-fw"></i>
<span>
Labels
</span>
</a></li><li class=""><a title="Wiki" class="shortcuts-wiki" href="/ftract_dev/ImaGIN2/wikis/home"><i class="fa fa-book fa-fw"></i>
<span>
Wiki
</span>
</a></li><li class=""><a title="Forks" href="/ftract_dev/ImaGIN2/forks"><i class="fa fa-code-fork fa-fw"></i>
<span>
Forks
</span>
</a></li><li class="separate-item"><a title="Settings" href="/ftract_dev/ImaGIN2/edit"><i class="fa fa-cogs fa-fw"></i>
<span>
Settings
</span>
</a></li><li class='hidden'>
<a title="Network" class="shortcuts-network" href="/ftract_dev/ImaGIN2/network/c4081fa6f23394f3f85ab37ca70e4f12c0752677">Network
</a></li>
<li class='hidden'>
<a class="shortcuts-new-issue" href="/ftract_dev/ImaGIN2/issues/new">Create a new issue
</a></li>
</ul>

<div class='collapse-nav'>
<a class="toggle-nav-collapse" title="Open/Close" href="#"><i class="fa fa-angle-left"></i></a>

</div>
<a class="sidebar-user" title="Profile" href="/u/anthonyboyer"><img alt="Profile" class="avatar avatar s36" src="/assets/no_avatar-07eeb128b993e74003e8664cff0b8e1e7234cec0443766a6763df32ca3472c23.png" />
<div class='username'>
anthonyboyer
</div>
</a></div>
<div class='content-wrapper controls-dropdown-visible'>


<div class='flash-container'>
</div>


<div class='container-fluid container-limited'>
<div class='content'>
<div class='clearfix'>


<div class='tree-holder' id='tree-holder'>
<div class='nav-block'>
<div class='tree-ref-holder'>
<form class="project-refs-form" action="/ftract_dev/ImaGIN2/refs/switch" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
<select name="ref" id="ref" class="project-refs-select select2 select2-sm"><optgroup label="Branches"><option value="master">master</option></optgroup><optgroup label="Tags"></optgroup><optgroup label="Commit"><option selected="selected" value="c4081fa6f23394f3f85ab37ca70e4f12c0752677">c4081fa6f23394f3f85ab37ca70e4f12c0752677</option></optgroup></select>
<input type="hidden" name="destination" id="destination" value="blob" />
<input type="hidden" name="path" id="path" value="toolbox/ImaGIN_Electrode.m" />
</form>


</div>
<ul class='breadcrumb repo-breadcrumb'>
<li>
<a href="/ftract_dev/ImaGIN2/tree/c4081fa6f23394f3f85ab37ca70e4f12c0752677">ImaGIN2
</a></li>
<li>
<a href="/ftract_dev/ImaGIN2/tree/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox">toolbox</a>
</li>
<li>
<a href="/ftract_dev/ImaGIN2/blob/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m"><strong>
ImaGIN_Electrode.m
</strong>
</a></li>
</ul>
</div>
<ul class='blob-commit-info hidden-xs'>
<li class='commit js-toggle-container' id='commit-aedb8927'>
<div class='commit-row-title'>
<span class='item-title'>
<a class="commit-row-message" href="/ftract_dev/ImaGIN2/commit/aedb892754df912ab2d71946e9611e3f5ccdff85">Fix how to manage ', p, - in electrode label.</a>
</span>
<div class='pull-right'>
<button class="btn btn-clipboard" data-clipboard-text="aedb892754df912ab2d71946e9611e3f5ccdff85" type="button"><i class="fa fa-clipboard"></i></button>
<a class="commit_short_id" href="/ftract_dev/ImaGIN2/commit/aedb892754df912ab2d71946e9611e3f5ccdff85">aedb8927</a>
</div>
</div>
<div class='commit-row-info'>
by
<a class="commit-author-link has-tooltip" title="anthony.boyer.gin@univ-grenoble-alpes.fr" href="/u/anthonyboyer"><img class="avatar s24" width="24" alt="" src="/assets/no_avatar-07eeb128b993e74003e8664cff0b8e1e7234cec0443766a6763df32ca3472c23.png" /> <span class="commit-author-name">Anthony Boyer</span></a>
<div class='committed_ago'>
<time class="time_ago js-timeago js-timeago-pending" datetime="2021-02-23T14:19:38Z" title="Feb 23, 2021 2:19pm" data-toggle="tooltip" data-placement="top" data-container="body">2021-02-23 15:19:38 +0100</time><script>
//<![CDATA[
$('.js-timeago-pending').removeClass('js-timeago-pending').timeago()
//]]>
</script> &nbsp;
</div>
<a class="pull-right" href="/ftract_dev/ImaGIN2/tree/aedb892754df912ab2d71946e9611e3f5ccdff85">Browse Files</a>
</div>
</li>

</ul>
<div class='blob-content-holder' id='blob-content-holder'>
<article class='file-holder'>
<div class='file-title'>
<i class="fa fa-file-text-o fa-fw"></i>
<strong>
ImaGIN_Electrode.m
</strong>
<small>
17.2 KB
</small>
<div class='file-actions hidden-xs'>
<div class='btn-group tree-btn-group'>
<a class="btn btn-sm" target="_blank" href="/ftract_dev/ImaGIN2/raw/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m">Raw</a>
<a class="btn btn-sm" href="/ftract_dev/ImaGIN2/blame/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m">Blame</a>
<a class="btn btn-sm" href="/ftract_dev/ImaGIN2/commits/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m">History</a>
<a class="btn btn-sm" href="/ftract_dev/ImaGIN2/blob/c4081fa6f23394f3f85ab37ca70e4f12c0752677/toolbox/ImaGIN_Electrode.m">Permalink</a>
</div>
<div class='btn-group' role='group'>
<button name="button" type="submit" class="btn disabled has-tooltip btn-file-option" title="You can only edit files when you are on a branch" data-container="body">Edit</button>
<button name="button" type="submit" class="btn btn-default disabled has-tooltip" title="You can only replace files when you are on a branch" data-container="body">Replace</button>
<button name="button" type="submit" class="btn btn-remove disabled has-tooltip" title="You can only delete files when you are on a branch" data-container="body">Delete</button>
</div>

</div>
</div>
<div class='file-content code js-syntax-highlight'>
<div class='line-numbers'>
<a class='diff-line-num' data-line-number='1' href='#L1' id='L1'>
<i class="fa fa-link"></i>
1
</a>
<a class='diff-line-num' data-line-number='2' href='#L2' id='L2'>
<i class="fa fa-link"></i>
2
</a>
<a class='diff-line-num' data-line-number='3' href='#L3' id='L3'>
<i class="fa fa-link"></i>
3
</a>
<a class='diff-line-num' data-line-number='4' href='#L4' id='L4'>
<i class="fa fa-link"></i>
4
</a>
<a class='diff-line-num' data-line-number='5' href='#L5' id='L5'>
<i class="fa fa-link"></i>
5
</a>
<a class='diff-line-num' data-line-number='6' href='#L6' id='L6'>
<i class="fa fa-link"></i>
6
</a>
<a class='diff-line-num' data-line-number='7' href='#L7' id='L7'>
<i class="fa fa-link"></i>
7
</a>
<a class='diff-line-num' data-line-number='8' href='#L8' id='L8'>
<i class="fa fa-link"></i>
8
</a>
<a class='diff-line-num' data-line-number='9' href='#L9' id='L9'>
<i class="fa fa-link"></i>
9
</a>
<a class='diff-line-num' data-line-number='10' href='#L10' id='L10'>
<i class="fa fa-link"></i>
10
</a>
<a class='diff-line-num' data-line-number='11' href='#L11' id='L11'>
<i class="fa fa-link"></i>
11
</a>
<a class='diff-line-num' data-line-number='12' href='#L12' id='L12'>
<i class="fa fa-link"></i>
12
</a>
<a class='diff-line-num' data-line-number='13' href='#L13' id='L13'>
<i class="fa fa-link"></i>
13
</a>
<a class='diff-line-num' data-line-number='14' href='#L14' id='L14'>
<i class="fa fa-link"></i>
14
</a>
<a class='diff-line-num' data-line-number='15' href='#L15' id='L15'>
<i class="fa fa-link"></i>
15
</a>
<a class='diff-line-num' data-line-number='16' href='#L16' id='L16'>
<i class="fa fa-link"></i>
16
</a>
<a class='diff-line-num' data-line-number='17' href='#L17' id='L17'>
<i class="fa fa-link"></i>
17
</a>
<a class='diff-line-num' data-line-number='18' href='#L18' id='L18'>
<i class="fa fa-link"></i>
18
</a>
<a class='diff-line-num' data-line-number='19' href='#L19' id='L19'>
<i class="fa fa-link"></i>
19
</a>
<a class='diff-line-num' data-line-number='20' href='#L20' id='L20'>
<i class="fa fa-link"></i>
20
</a>
<a class='diff-line-num' data-line-number='21' href='#L21' id='L21'>
<i class="fa fa-link"></i>
21
</a>
<a class='diff-line-num' data-line-number='22' href='#L22' id='L22'>
<i class="fa fa-link"></i>
22
</a>
<a class='diff-line-num' data-line-number='23' href='#L23' id='L23'>
<i class="fa fa-link"></i>
23
</a>
<a class='diff-line-num' data-line-number='24' href='#L24' id='L24'>
<i class="fa fa-link"></i>
24
</a>
<a class='diff-line-num' data-line-number='25' href='#L25' id='L25'>
<i class="fa fa-link"></i>
25
</a>
<a class='diff-line-num' data-line-number='26' href='#L26' id='L26'>
<i class="fa fa-link"></i>
26
</a>
<a class='diff-line-num' data-line-number='27' href='#L27' id='L27'>
<i class="fa fa-link"></i>
27
</a>
<a class='diff-line-num' data-line-number='28' href='#L28' id='L28'>
<i class="fa fa-link"></i>
28
</a>
<a class='diff-line-num' data-line-number='29' href='#L29' id='L29'>
<i class="fa fa-link"></i>
29
</a>
<a class='diff-line-num' data-line-number='30' href='#L30' id='L30'>
<i class="fa fa-link"></i>
30
</a>
<a class='diff-line-num' data-line-number='31' href='#L31' id='L31'>
<i class="fa fa-link"></i>
31
</a>
<a class='diff-line-num' data-line-number='32' href='#L32' id='L32'>
<i class="fa fa-link"></i>
32
</a>
<a class='diff-line-num' data-line-number='33' href='#L33' id='L33'>
<i class="fa fa-link"></i>
33
</a>
<a class='diff-line-num' data-line-number='34' href='#L34' id='L34'>
<i class="fa fa-link"></i>
34
</a>
<a class='diff-line-num' data-line-number='35' href='#L35' id='L35'>
<i class="fa fa-link"></i>
35
</a>
<a class='diff-line-num' data-line-number='36' href='#L36' id='L36'>
<i class="fa fa-link"></i>
36
</a>
<a class='diff-line-num' data-line-number='37' href='#L37' id='L37'>
<i class="fa fa-link"></i>
37
</a>
<a class='diff-line-num' data-line-number='38' href='#L38' id='L38'>
<i class="fa fa-link"></i>
38
</a>
<a class='diff-line-num' data-line-number='39' href='#L39' id='L39'>
<i class="fa fa-link"></i>
39
</a>
<a class='diff-line-num' data-line-number='40' href='#L40' id='L40'>
<i class="fa fa-link"></i>
40
</a>
<a class='diff-line-num' data-line-number='41' href='#L41' id='L41'>
<i class="fa fa-link"></i>
41
</a>
<a class='diff-line-num' data-line-number='42' href='#L42' id='L42'>
<i class="fa fa-link"></i>
42
</a>
<a class='diff-line-num' data-line-number='43' href='#L43' id='L43'>
<i class="fa fa-link"></i>
43
</a>
<a class='diff-line-num' data-line-number='44' href='#L44' id='L44'>
<i class="fa fa-link"></i>
44
</a>
<a class='diff-line-num' data-line-number='45' href='#L45' id='L45'>
<i class="fa fa-link"></i>
45
</a>
<a class='diff-line-num' data-line-number='46' href='#L46' id='L46'>
<i class="fa fa-link"></i>
46
</a>
<a class='diff-line-num' data-line-number='47' href='#L47' id='L47'>
<i class="fa fa-link"></i>
47
</a>
<a class='diff-line-num' data-line-number='48' href='#L48' id='L48'>
<i class="fa fa-link"></i>
48
</a>
<a class='diff-line-num' data-line-number='49' href='#L49' id='L49'>
<i class="fa fa-link"></i>
49
</a>
<a class='diff-line-num' data-line-number='50' href='#L50' id='L50'>
<i class="fa fa-link"></i>
50
</a>
<a class='diff-line-num' data-line-number='51' href='#L51' id='L51'>
<i class="fa fa-link"></i>
51
</a>
<a class='diff-line-num' data-line-number='52' href='#L52' id='L52'>
<i class="fa fa-link"></i>
52
</a>
<a class='diff-line-num' data-line-number='53' href='#L53' id='L53'>
<i class="fa fa-link"></i>
53
</a>
<a class='diff-line-num' data-line-number='54' href='#L54' id='L54'>
<i class="fa fa-link"></i>
54
</a>
<a class='diff-line-num' data-line-number='55' href='#L55' id='L55'>
<i class="fa fa-link"></i>
55
</a>
<a class='diff-line-num' data-line-number='56' href='#L56' id='L56'>
<i class="fa fa-link"></i>
56
</a>
<a class='diff-line-num' data-line-number='57' href='#L57' id='L57'>
<i class="fa fa-link"></i>
57
</a>
<a class='diff-line-num' data-line-number='58' href='#L58' id='L58'>
<i class="fa fa-link"></i>
58
</a>
<a class='diff-line-num' data-line-number='59' href='#L59' id='L59'>
<i class="fa fa-link"></i>
59
</a>
<a class='diff-line-num' data-line-number='60' href='#L60' id='L60'>
<i class="fa fa-link"></i>
60
</a>
<a class='diff-line-num' data-line-number='61' href='#L61' id='L61'>
<i class="fa fa-link"></i>
61
</a>
<a class='diff-line-num' data-line-number='62' href='#L62' id='L62'>
<i class="fa fa-link"></i>
62
</a>
<a class='diff-line-num' data-line-number='63' href='#L63' id='L63'>
<i class="fa fa-link"></i>
63
</a>
<a class='diff-line-num' data-line-number='64' href='#L64' id='L64'>
<i class="fa fa-link"></i>
64
</a>
<a class='diff-line-num' data-line-number='65' href='#L65' id='L65'>
<i class="fa fa-link"></i>
65
</a>
<a class='diff-line-num' data-line-number='66' href='#L66' id='L66'>
<i class="fa fa-link"></i>
66
</a>
<a class='diff-line-num' data-line-number='67' href='#L67' id='L67'>
<i class="fa fa-link"></i>
67
</a>
<a class='diff-line-num' data-line-number='68' href='#L68' id='L68'>
<i class="fa fa-link"></i>
68
</a>
<a class='diff-line-num' data-line-number='69' href='#L69' id='L69'>
<i class="fa fa-link"></i>
69
</a>
<a class='diff-line-num' data-line-number='70' href='#L70' id='L70'>
<i class="fa fa-link"></i>
70
</a>
<a class='diff-line-num' data-line-number='71' href='#L71' id='L71'>
<i class="fa fa-link"></i>
71
</a>
<a class='diff-line-num' data-line-number='72' href='#L72' id='L72'>
<i class="fa fa-link"></i>
72
</a>
<a class='diff-line-num' data-line-number='73' href='#L73' id='L73'>
<i class="fa fa-link"></i>
73
</a>
<a class='diff-line-num' data-line-number='74' href='#L74' id='L74'>
<i class="fa fa-link"></i>
74
</a>
<a class='diff-line-num' data-line-number='75' href='#L75' id='L75'>
<i class="fa fa-link"></i>
75
</a>
<a class='diff-line-num' data-line-number='76' href='#L76' id='L76'>
<i class="fa fa-link"></i>
76
</a>
<a class='diff-line-num' data-line-number='77' href='#L77' id='L77'>
<i class="fa fa-link"></i>
77
</a>
<a class='diff-line-num' data-line-number='78' href='#L78' id='L78'>
<i class="fa fa-link"></i>
78
</a>
<a class='diff-line-num' data-line-number='79' href='#L79' id='L79'>
<i class="fa fa-link"></i>
79
</a>
<a class='diff-line-num' data-line-number='80' href='#L80' id='L80'>
<i class="fa fa-link"></i>
80
</a>
<a class='diff-line-num' data-line-number='81' href='#L81' id='L81'>
<i class="fa fa-link"></i>
81
</a>
<a class='diff-line-num' data-line-number='82' href='#L82' id='L82'>
<i class="fa fa-link"></i>
82
</a>
<a class='diff-line-num' data-line-number='83' href='#L83' id='L83'>
<i class="fa fa-link"></i>
83
</a>
<a class='diff-line-num' data-line-number='84' href='#L84' id='L84'>
<i class="fa fa-link"></i>
84
</a>
<a class='diff-line-num' data-line-number='85' href='#L85' id='L85'>
<i class="fa fa-link"></i>
85
</a>
<a class='diff-line-num' data-line-number='86' href='#L86' id='L86'>
<i class="fa fa-link"></i>
86
</a>
<a class='diff-line-num' data-line-number='87' href='#L87' id='L87'>
<i class="fa fa-link"></i>
87
</a>
<a class='diff-line-num' data-line-number='88' href='#L88' id='L88'>
<i class="fa fa-link"></i>
88
</a>
<a class='diff-line-num' data-line-number='89' href='#L89' id='L89'>
<i class="fa fa-link"></i>
89
</a>
<a class='diff-line-num' data-line-number='90' href='#L90' id='L90'>
<i class="fa fa-link"></i>
90
</a>
<a class='diff-line-num' data-line-number='91' href='#L91' id='L91'>
<i class="fa fa-link"></i>
91
</a>
<a class='diff-line-num' data-line-number='92' href='#L92' id='L92'>
<i class="fa fa-link"></i>
92
</a>
<a class='diff-line-num' data-line-number='93' href='#L93' id='L93'>
<i class="fa fa-link"></i>
93
</a>
<a class='diff-line-num' data-line-number='94' href='#L94' id='L94'>
<i class="fa fa-link"></i>
94
</a>
<a class='diff-line-num' data-line-number='95' href='#L95' id='L95'>
<i class="fa fa-link"></i>
95
</a>
<a class='diff-line-num' data-line-number='96' href='#L96' id='L96'>
<i class="fa fa-link"></i>
96
</a>
<a class='diff-line-num' data-line-number='97' href='#L97' id='L97'>
<i class="fa fa-link"></i>
97
</a>
<a class='diff-line-num' data-line-number='98' href='#L98' id='L98'>
<i class="fa fa-link"></i>
98
</a>
<a class='diff-line-num' data-line-number='99' href='#L99' id='L99'>
<i class="fa fa-link"></i>
99
</a>
<a class='diff-line-num' data-line-number='100' href='#L100' id='L100'>
<i class="fa fa-link"></i>
100
</a>
<a class='diff-line-num' data-line-number='101' href='#L101' id='L101'>
<i class="fa fa-link"></i>
101
</a>
<a class='diff-line-num' data-line-number='102' href='#L102' id='L102'>
<i class="fa fa-link"></i>
102
</a>
<a class='diff-line-num' data-line-number='103' href='#L103' id='L103'>
<i class="fa fa-link"></i>
103
</a>
<a class='diff-line-num' data-line-number='104' href='#L104' id='L104'>
<i class="fa fa-link"></i>
104
</a>
<a class='diff-line-num' data-line-number='105' href='#L105' id='L105'>
<i class="fa fa-link"></i>
105
</a>
<a class='diff-line-num' data-line-number='106' href='#L106' id='L106'>
<i class="fa fa-link"></i>
106
</a>
<a class='diff-line-num' data-line-number='107' href='#L107' id='L107'>
<i class="fa fa-link"></i>
107
</a>
<a class='diff-line-num' data-line-number='108' href='#L108' id='L108'>
<i class="fa fa-link"></i>
108
</a>
<a class='diff-line-num' data-line-number='109' href='#L109' id='L109'>
<i class="fa fa-link"></i>
109
</a>
<a class='diff-line-num' data-line-number='110' href='#L110' id='L110'>
<i class="fa fa-link"></i>
110
</a>
<a class='diff-line-num' data-line-number='111' href='#L111' id='L111'>
<i class="fa fa-link"></i>
111
</a>
<a class='diff-line-num' data-line-number='112' href='#L112' id='L112'>
<i class="fa fa-link"></i>
112
</a>
<a class='diff-line-num' data-line-number='113' href='#L113' id='L113'>
<i class="fa fa-link"></i>
113
</a>
<a class='diff-line-num' data-line-number='114' href='#L114' id='L114'>
<i class="fa fa-link"></i>
114
</a>
<a class='diff-line-num' data-line-number='115' href='#L115' id='L115'>
<i class="fa fa-link"></i>
115
</a>
<a class='diff-line-num' data-line-number='116' href='#L116' id='L116'>
<i class="fa fa-link"></i>
116
</a>
<a class='diff-line-num' data-line-number='117' href='#L117' id='L117'>
<i class="fa fa-link"></i>
117
</a>
<a class='diff-line-num' data-line-number='118' href='#L118' id='L118'>
<i class="fa fa-link"></i>
118
</a>
<a class='diff-line-num' data-line-number='119' href='#L119' id='L119'>
<i class="fa fa-link"></i>
119
</a>
<a class='diff-line-num' data-line-number='120' href='#L120' id='L120'>
<i class="fa fa-link"></i>
120
</a>
<a class='diff-line-num' data-line-number='121' href='#L121' id='L121'>
<i class="fa fa-link"></i>
121
</a>
<a class='diff-line-num' data-line-number='122' href='#L122' id='L122'>
<i class="fa fa-link"></i>
122
</a>
<a class='diff-line-num' data-line-number='123' href='#L123' id='L123'>
<i class="fa fa-link"></i>
123
</a>
<a class='diff-line-num' data-line-number='124' href='#L124' id='L124'>
<i class="fa fa-link"></i>
124
</a>
<a class='diff-line-num' data-line-number='125' href='#L125' id='L125'>
<i class="fa fa-link"></i>
125
</a>
<a class='diff-line-num' data-line-number='126' href='#L126' id='L126'>
<i class="fa fa-link"></i>
126
</a>
<a class='diff-line-num' data-line-number='127' href='#L127' id='L127'>
<i class="fa fa-link"></i>
127
</a>
<a class='diff-line-num' data-line-number='128' href='#L128' id='L128'>
<i class="fa fa-link"></i>
128
</a>
<a class='diff-line-num' data-line-number='129' href='#L129' id='L129'>
<i class="fa fa-link"></i>
129
</a>
<a class='diff-line-num' data-line-number='130' href='#L130' id='L130'>
<i class="fa fa-link"></i>
130
</a>
<a class='diff-line-num' data-line-number='131' href='#L131' id='L131'>
<i class="fa fa-link"></i>
131
</a>
<a class='diff-line-num' data-line-number='132' href='#L132' id='L132'>
<i class="fa fa-link"></i>
132
</a>
<a class='diff-line-num' data-line-number='133' href='#L133' id='L133'>
<i class="fa fa-link"></i>
133
</a>
<a class='diff-line-num' data-line-number='134' href='#L134' id='L134'>
<i class="fa fa-link"></i>
134
</a>
<a class='diff-line-num' data-line-number='135' href='#L135' id='L135'>
<i class="fa fa-link"></i>
135
</a>
<a class='diff-line-num' data-line-number='136' href='#L136' id='L136'>
<i class="fa fa-link"></i>
136
</a>
<a class='diff-line-num' data-line-number='137' href='#L137' id='L137'>
<i class="fa fa-link"></i>
137
</a>
<a class='diff-line-num' data-line-number='138' href='#L138' id='L138'>
<i class="fa fa-link"></i>
138
</a>
<a class='diff-line-num' data-line-number='139' href='#L139' id='L139'>
<i class="fa fa-link"></i>
139
</a>
<a class='diff-line-num' data-line-number='140' href='#L140' id='L140'>
<i class="fa fa-link"></i>
140
</a>
<a class='diff-line-num' data-line-number='141' href='#L141' id='L141'>
<i class="fa fa-link"></i>
141
</a>
<a class='diff-line-num' data-line-number='142' href='#L142' id='L142'>
<i class="fa fa-link"></i>
142
</a>
<a class='diff-line-num' data-line-number='143' href='#L143' id='L143'>
<i class="fa fa-link"></i>
143
</a>
<a class='diff-line-num' data-line-number='144' href='#L144' id='L144'>
<i class="fa fa-link"></i>
144
</a>
<a class='diff-line-num' data-line-number='145' href='#L145' id='L145'>
<i class="fa fa-link"></i>
145
</a>
<a class='diff-line-num' data-line-number='146' href='#L146' id='L146'>
<i class="fa fa-link"></i>
146
</a>
<a class='diff-line-num' data-line-number='147' href='#L147' id='L147'>
<i class="fa fa-link"></i>
147
</a>
<a class='diff-line-num' data-line-number='148' href='#L148' id='L148'>
<i class="fa fa-link"></i>
148
</a>
<a class='diff-line-num' data-line-number='149' href='#L149' id='L149'>
<i class="fa fa-link"></i>
149
</a>
<a class='diff-line-num' data-line-number='150' href='#L150' id='L150'>
<i class="fa fa-link"></i>
150
</a>
<a class='diff-line-num' data-line-number='151' href='#L151' id='L151'>
<i class="fa fa-link"></i>
151
</a>
<a class='diff-line-num' data-line-number='152' href='#L152' id='L152'>
<i class="fa fa-link"></i>
152
</a>
<a class='diff-line-num' data-line-number='153' href='#L153' id='L153'>
<i class="fa fa-link"></i>
153
</a>
<a class='diff-line-num' data-line-number='154' href='#L154' id='L154'>
<i class="fa fa-link"></i>
154
</a>
<a class='diff-line-num' data-line-number='155' href='#L155' id='L155'>
<i class="fa fa-link"></i>
155
</a>
<a class='diff-line-num' data-line-number='156' href='#L156' id='L156'>
<i class="fa fa-link"></i>
156
</a>
<a class='diff-line-num' data-line-number='157' href='#L157' id='L157'>
<i class="fa fa-link"></i>
157
</a>
<a class='diff-line-num' data-line-number='158' href='#L158' id='L158'>
<i class="fa fa-link"></i>
158
</a>
<a class='diff-line-num' data-line-number='159' href='#L159' id='L159'>
<i class="fa fa-link"></i>
159
</a>
<a class='diff-line-num' data-line-number='160' href='#L160' id='L160'>
<i class="fa fa-link"></i>
160
</a>
<a class='diff-line-num' data-line-number='161' href='#L161' id='L161'>
<i class="fa fa-link"></i>
161
</a>
<a class='diff-line-num' data-line-number='162' href='#L162' id='L162'>
<i class="fa fa-link"></i>
162
</a>
<a class='diff-line-num' data-line-number='163' href='#L163' id='L163'>
<i class="fa fa-link"></i>
163
</a>
<a class='diff-line-num' data-line-number='164' href='#L164' id='L164'>
<i class="fa fa-link"></i>
164
</a>
<a class='diff-line-num' data-line-number='165' href='#L165' id='L165'>
<i class="fa fa-link"></i>
165
</a>
<a class='diff-line-num' data-line-number='166' href='#L166' id='L166'>
<i class="fa fa-link"></i>
166
</a>
<a class='diff-line-num' data-line-number='167' href='#L167' id='L167'>
<i class="fa fa-link"></i>
167
</a>
<a class='diff-line-num' data-line-number='168' href='#L168' id='L168'>
<i class="fa fa-link"></i>
168
</a>
<a class='diff-line-num' data-line-number='169' href='#L169' id='L169'>
<i class="fa fa-link"></i>
169
</a>
<a class='diff-line-num' data-line-number='170' href='#L170' id='L170'>
<i class="fa fa-link"></i>
170
</a>
<a class='diff-line-num' data-line-number='171' href='#L171' id='L171'>
<i class="fa fa-link"></i>
171
</a>
<a class='diff-line-num' data-line-number='172' href='#L172' id='L172'>
<i class="fa fa-link"></i>
172
</a>
<a class='diff-line-num' data-line-number='173' href='#L173' id='L173'>
<i class="fa fa-link"></i>
173
</a>
<a class='diff-line-num' data-line-number='174' href='#L174' id='L174'>
<i class="fa fa-link"></i>
174
</a>
<a class='diff-line-num' data-line-number='175' href='#L175' id='L175'>
<i class="fa fa-link"></i>
175
</a>
<a class='diff-line-num' data-line-number='176' href='#L176' id='L176'>
<i class="fa fa-link"></i>
176
</a>
<a class='diff-line-num' data-line-number='177' href='#L177' id='L177'>
<i class="fa fa-link"></i>
177
</a>
<a class='diff-line-num' data-line-number='178' href='#L178' id='L178'>
<i class="fa fa-link"></i>
178
</a>
<a class='diff-line-num' data-line-number='179' href='#L179' id='L179'>
<i class="fa fa-link"></i>
179
</a>
<a class='diff-line-num' data-line-number='180' href='#L180' id='L180'>
<i class="fa fa-link"></i>
180
</a>
<a class='diff-line-num' data-line-number='181' href='#L181' id='L181'>
<i class="fa fa-link"></i>
181
</a>
<a class='diff-line-num' data-line-number='182' href='#L182' id='L182'>
<i class="fa fa-link"></i>
182
</a>
<a class='diff-line-num' data-line-number='183' href='#L183' id='L183'>
<i class="fa fa-link"></i>
183
</a>
<a class='diff-line-num' data-line-number='184' href='#L184' id='L184'>
<i class="fa fa-link"></i>
184
</a>
<a class='diff-line-num' data-line-number='185' href='#L185' id='L185'>
<i class="fa fa-link"></i>
185
</a>
<a class='diff-line-num' data-line-number='186' href='#L186' id='L186'>
<i class="fa fa-link"></i>
186
</a>
<a class='diff-line-num' data-line-number='187' href='#L187' id='L187'>
<i class="fa fa-link"></i>
187
</a>
<a class='diff-line-num' data-line-number='188' href='#L188' id='L188'>
<i class="fa fa-link"></i>
188
</a>
<a class='diff-line-num' data-line-number='189' href='#L189' id='L189'>
<i class="fa fa-link"></i>
189
</a>
<a class='diff-line-num' data-line-number='190' href='#L190' id='L190'>
<i class="fa fa-link"></i>
190
</a>
<a class='diff-line-num' data-line-number='191' href='#L191' id='L191'>
<i class="fa fa-link"></i>
191
</a>
<a class='diff-line-num' data-line-number='192' href='#L192' id='L192'>
<i class="fa fa-link"></i>
192
</a>
<a class='diff-line-num' data-line-number='193' href='#L193' id='L193'>
<i class="fa fa-link"></i>
193
</a>
<a class='diff-line-num' data-line-number='194' href='#L194' id='L194'>
<i class="fa fa-link"></i>
194
</a>
<a class='diff-line-num' data-line-number='195' href='#L195' id='L195'>
<i class="fa fa-link"></i>
195
</a>
<a class='diff-line-num' data-line-number='196' href='#L196' id='L196'>
<i class="fa fa-link"></i>
196
</a>
<a class='diff-line-num' data-line-number='197' href='#L197' id='L197'>
<i class="fa fa-link"></i>
197
</a>
<a class='diff-line-num' data-line-number='198' href='#L198' id='L198'>
<i class="fa fa-link"></i>
198
</a>
<a class='diff-line-num' data-line-number='199' href='#L199' id='L199'>
<i class="fa fa-link"></i>
199
</a>
<a class='diff-line-num' data-line-number='200' href='#L200' id='L200'>
<i class="fa fa-link"></i>
200
</a>
<a class='diff-line-num' data-line-number='201' href='#L201' id='L201'>
<i class="fa fa-link"></i>
201
</a>
<a class='diff-line-num' data-line-number='202' href='#L202' id='L202'>
<i class="fa fa-link"></i>
202
</a>
<a class='diff-line-num' data-line-number='203' href='#L203' id='L203'>
<i class="fa fa-link"></i>
203
</a>
<a class='diff-line-num' data-line-number='204' href='#L204' id='L204'>
<i class="fa fa-link"></i>
204
</a>
<a class='diff-line-num' data-line-number='205' href='#L205' id='L205'>
<i class="fa fa-link"></i>
205
</a>
<a class='diff-line-num' data-line-number='206' href='#L206' id='L206'>
<i class="fa fa-link"></i>
206
</a>
<a class='diff-line-num' data-line-number='207' href='#L207' id='L207'>
<i class="fa fa-link"></i>
207
</a>
<a class='diff-line-num' data-line-number='208' href='#L208' id='L208'>
<i class="fa fa-link"></i>
208
</a>
<a class='diff-line-num' data-line-number='209' href='#L209' id='L209'>
<i class="fa fa-link"></i>
209
</a>
<a class='diff-line-num' data-line-number='210' href='#L210' id='L210'>
<i class="fa fa-link"></i>
210
</a>
<a class='diff-line-num' data-line-number='211' href='#L211' id='L211'>
<i class="fa fa-link"></i>
211
</a>
<a class='diff-line-num' data-line-number='212' href='#L212' id='L212'>
<i class="fa fa-link"></i>
212
</a>
<a class='diff-line-num' data-line-number='213' href='#L213' id='L213'>
<i class="fa fa-link"></i>
213
</a>
<a class='diff-line-num' data-line-number='214' href='#L214' id='L214'>
<i class="fa fa-link"></i>
214
</a>
<a class='diff-line-num' data-line-number='215' href='#L215' id='L215'>
<i class="fa fa-link"></i>
215
</a>
<a class='diff-line-num' data-line-number='216' href='#L216' id='L216'>
<i class="fa fa-link"></i>
216
</a>
<a class='diff-line-num' data-line-number='217' href='#L217' id='L217'>
<i class="fa fa-link"></i>
217
</a>
<a class='diff-line-num' data-line-number='218' href='#L218' id='L218'>
<i class="fa fa-link"></i>
218
</a>
<a class='diff-line-num' data-line-number='219' href='#L219' id='L219'>
<i class="fa fa-link"></i>
219
</a>
<a class='diff-line-num' data-line-number='220' href='#L220' id='L220'>
<i class="fa fa-link"></i>
220
</a>
<a class='diff-line-num' data-line-number='221' href='#L221' id='L221'>
<i class="fa fa-link"></i>
221
</a>
<a class='diff-line-num' data-line-number='222' href='#L222' id='L222'>
<i class="fa fa-link"></i>
222
</a>
<a class='diff-line-num' data-line-number='223' href='#L223' id='L223'>
<i class="fa fa-link"></i>
223
</a>
<a class='diff-line-num' data-line-number='224' href='#L224' id='L224'>
<i class="fa fa-link"></i>
224
</a>
<a class='diff-line-num' data-line-number='225' href='#L225' id='L225'>
<i class="fa fa-link"></i>
225
</a>
<a class='diff-line-num' data-line-number='226' href='#L226' id='L226'>
<i class="fa fa-link"></i>
226
</a>
<a class='diff-line-num' data-line-number='227' href='#L227' id='L227'>
<i class="fa fa-link"></i>
227
</a>
<a class='diff-line-num' data-line-number='228' href='#L228' id='L228'>
<i class="fa fa-link"></i>
228
</a>
<a class='diff-line-num' data-line-number='229' href='#L229' id='L229'>
<i class="fa fa-link"></i>
229
</a>
<a class='diff-line-num' data-line-number='230' href='#L230' id='L230'>
<i class="fa fa-link"></i>
230
</a>
<a class='diff-line-num' data-line-number='231' href='#L231' id='L231'>
<i class="fa fa-link"></i>
231
</a>
<a class='diff-line-num' data-line-number='232' href='#L232' id='L232'>
<i class="fa fa-link"></i>
232
</a>
<a class='diff-line-num' data-line-number='233' href='#L233' id='L233'>
<i class="fa fa-link"></i>
233
</a>
<a class='diff-line-num' data-line-number='234' href='#L234' id='L234'>
<i class="fa fa-link"></i>
234
</a>
<a class='diff-line-num' data-line-number='235' href='#L235' id='L235'>
<i class="fa fa-link"></i>
235
</a>
<a class='diff-line-num' data-line-number='236' href='#L236' id='L236'>
<i class="fa fa-link"></i>
236
</a>
<a class='diff-line-num' data-line-number='237' href='#L237' id='L237'>
<i class="fa fa-link"></i>
237
</a>
<a class='diff-line-num' data-line-number='238' href='#L238' id='L238'>
<i class="fa fa-link"></i>
238
</a>
<a class='diff-line-num' data-line-number='239' href='#L239' id='L239'>
<i class="fa fa-link"></i>
239
</a>
<a class='diff-line-num' data-line-number='240' href='#L240' id='L240'>
<i class="fa fa-link"></i>
240
</a>
<a class='diff-line-num' data-line-number='241' href='#L241' id='L241'>
<i class="fa fa-link"></i>
241
</a>
<a class='diff-line-num' data-line-number='242' href='#L242' id='L242'>
<i class="fa fa-link"></i>
242
</a>
<a class='diff-line-num' data-line-number='243' href='#L243' id='L243'>
<i class="fa fa-link"></i>
243
</a>
<a class='diff-line-num' data-line-number='244' href='#L244' id='L244'>
<i class="fa fa-link"></i>
244
</a>
<a class='diff-line-num' data-line-number='245' href='#L245' id='L245'>
<i class="fa fa-link"></i>
245
</a>
<a class='diff-line-num' data-line-number='246' href='#L246' id='L246'>
<i class="fa fa-link"></i>
246
</a>
<a class='diff-line-num' data-line-number='247' href='#L247' id='L247'>
<i class="fa fa-link"></i>
247
</a>
<a class='diff-line-num' data-line-number='248' href='#L248' id='L248'>
<i class="fa fa-link"></i>
248
</a>
<a class='diff-line-num' data-line-number='249' href='#L249' id='L249'>
<i class="fa fa-link"></i>
249
</a>
<a class='diff-line-num' data-line-number='250' href='#L250' id='L250'>
<i class="fa fa-link"></i>
250
</a>
<a class='diff-line-num' data-line-number='251' href='#L251' id='L251'>
<i class="fa fa-link"></i>
251
</a>
<a class='diff-line-num' data-line-number='252' href='#L252' id='L252'>
<i class="fa fa-link"></i>
252
</a>
<a class='diff-line-num' data-line-number='253' href='#L253' id='L253'>
<i class="fa fa-link"></i>
253
</a>
<a class='diff-line-num' data-line-number='254' href='#L254' id='L254'>
<i class="fa fa-link"></i>
254
</a>
<a class='diff-line-num' data-line-number='255' href='#L255' id='L255'>
<i class="fa fa-link"></i>
255
</a>
<a class='diff-line-num' data-line-number='256' href='#L256' id='L256'>
<i class="fa fa-link"></i>
256
</a>
<a class='diff-line-num' data-line-number='257' href='#L257' id='L257'>
<i class="fa fa-link"></i>
257
</a>
<a class='diff-line-num' data-line-number='258' href='#L258' id='L258'>
<i class="fa fa-link"></i>
258
</a>
<a class='diff-line-num' data-line-number='259' href='#L259' id='L259'>
<i class="fa fa-link"></i>
259
</a>
<a class='diff-line-num' data-line-number='260' href='#L260' id='L260'>
<i class="fa fa-link"></i>
260
</a>
<a class='diff-line-num' data-line-number='261' href='#L261' id='L261'>
<i class="fa fa-link"></i>
261
</a>
<a class='diff-line-num' data-line-number='262' href='#L262' id='L262'>
<i class="fa fa-link"></i>
262
</a>
<a class='diff-line-num' data-line-number='263' href='#L263' id='L263'>
<i class="fa fa-link"></i>
263
</a>
<a class='diff-line-num' data-line-number='264' href='#L264' id='L264'>
<i class="fa fa-link"></i>
264
</a>
<a class='diff-line-num' data-line-number='265' href='#L265' id='L265'>
<i class="fa fa-link"></i>
265
</a>
<a class='diff-line-num' data-line-number='266' href='#L266' id='L266'>
<i class="fa fa-link"></i>
266
</a>
<a class='diff-line-num' data-line-number='267' href='#L267' id='L267'>
<i class="fa fa-link"></i>
267
</a>
<a class='diff-line-num' data-line-number='268' href='#L268' id='L268'>
<i class="fa fa-link"></i>
268
</a>
<a class='diff-line-num' data-line-number='269' href='#L269' id='L269'>
<i class="fa fa-link"></i>
269
</a>
<a class='diff-line-num' data-line-number='270' href='#L270' id='L270'>
<i class="fa fa-link"></i>
270
</a>
<a class='diff-line-num' data-line-number='271' href='#L271' id='L271'>
<i class="fa fa-link"></i>
271
</a>
<a class='diff-line-num' data-line-number='272' href='#L272' id='L272'>
<i class="fa fa-link"></i>
272
</a>
<a class='diff-line-num' data-line-number='273' href='#L273' id='L273'>
<i class="fa fa-link"></i>
273
</a>
<a class='diff-line-num' data-line-number='274' href='#L274' id='L274'>
<i class="fa fa-link"></i>
274
</a>
<a class='diff-line-num' data-line-number='275' href='#L275' id='L275'>
<i class="fa fa-link"></i>
275
</a>
<a class='diff-line-num' data-line-number='276' href='#L276' id='L276'>
<i class="fa fa-link"></i>
276
</a>
<a class='diff-line-num' data-line-number='277' href='#L277' id='L277'>
<i class="fa fa-link"></i>
277
</a>
<a class='diff-line-num' data-line-number='278' href='#L278' id='L278'>
<i class="fa fa-link"></i>
278
</a>
<a class='diff-line-num' data-line-number='279' href='#L279' id='L279'>
<i class="fa fa-link"></i>
279
</a>
<a class='diff-line-num' data-line-number='280' href='#L280' id='L280'>
<i class="fa fa-link"></i>
280
</a>
<a class='diff-line-num' data-line-number='281' href='#L281' id='L281'>
<i class="fa fa-link"></i>
281
</a>
<a class='diff-line-num' data-line-number='282' href='#L282' id='L282'>
<i class="fa fa-link"></i>
282
</a>
<a class='diff-line-num' data-line-number='283' href='#L283' id='L283'>
<i class="fa fa-link"></i>
283
</a>
<a class='diff-line-num' data-line-number='284' href='#L284' id='L284'>
<i class="fa fa-link"></i>
284
</a>
<a class='diff-line-num' data-line-number='285' href='#L285' id='L285'>
<i class="fa fa-link"></i>
285
</a>
<a class='diff-line-num' data-line-number='286' href='#L286' id='L286'>
<i class="fa fa-link"></i>
286
</a>
<a class='diff-line-num' data-line-number='287' href='#L287' id='L287'>
<i class="fa fa-link"></i>
287
</a>
<a class='diff-line-num' data-line-number='288' href='#L288' id='L288'>
<i class="fa fa-link"></i>
288
</a>
<a class='diff-line-num' data-line-number='289' href='#L289' id='L289'>
<i class="fa fa-link"></i>
289
</a>
<a class='diff-line-num' data-line-number='290' href='#L290' id='L290'>
<i class="fa fa-link"></i>
290
</a>
<a class='diff-line-num' data-line-number='291' href='#L291' id='L291'>
<i class="fa fa-link"></i>
291
</a>
<a class='diff-line-num' data-line-number='292' href='#L292' id='L292'>
<i class="fa fa-link"></i>
292
</a>
<a class='diff-line-num' data-line-number='293' href='#L293' id='L293'>
<i class="fa fa-link"></i>
293
</a>
<a class='diff-line-num' data-line-number='294' href='#L294' id='L294'>
<i class="fa fa-link"></i>
294
</a>
<a class='diff-line-num' data-line-number='295' href='#L295' id='L295'>
<i class="fa fa-link"></i>
295
</a>
<a class='diff-line-num' data-line-number='296' href='#L296' id='L296'>
<i class="fa fa-link"></i>
296
</a>
<a class='diff-line-num' data-line-number='297' href='#L297' id='L297'>
<i class="fa fa-link"></i>
297
</a>
<a class='diff-line-num' data-line-number='298' href='#L298' id='L298'>
<i class="fa fa-link"></i>
298
</a>
<a class='diff-line-num' data-line-number='299' href='#L299' id='L299'>
<i class="fa fa-link"></i>
299
</a>
<a class='diff-line-num' data-line-number='300' href='#L300' id='L300'>
<i class="fa fa-link"></i>
300
</a>
<a class='diff-line-num' data-line-number='301' href='#L301' id='L301'>
<i class="fa fa-link"></i>
301
</a>
<a class='diff-line-num' data-line-number='302' href='#L302' id='L302'>
<i class="fa fa-link"></i>
302
</a>
<a class='diff-line-num' data-line-number='303' href='#L303' id='L303'>
<i class="fa fa-link"></i>
303
</a>
<a class='diff-line-num' data-line-number='304' href='#L304' id='L304'>
<i class="fa fa-link"></i>
304
</a>
<a class='diff-line-num' data-line-number='305' href='#L305' id='L305'>
<i class="fa fa-link"></i>
305
</a>
<a class='diff-line-num' data-line-number='306' href='#L306' id='L306'>
<i class="fa fa-link"></i>
306
</a>
<a class='diff-line-num' data-line-number='307' href='#L307' id='L307'>
<i class="fa fa-link"></i>
307
</a>
<a class='diff-line-num' data-line-number='308' href='#L308' id='L308'>
<i class="fa fa-link"></i>
308
</a>
<a class='diff-line-num' data-line-number='309' href='#L309' id='L309'>
<i class="fa fa-link"></i>
309
</a>
<a class='diff-line-num' data-line-number='310' href='#L310' id='L310'>
<i class="fa fa-link"></i>
310
</a>
<a class='diff-line-num' data-line-number='311' href='#L311' id='L311'>
<i class="fa fa-link"></i>
311
</a>
<a class='diff-line-num' data-line-number='312' href='#L312' id='L312'>
<i class="fa fa-link"></i>
312
</a>
<a class='diff-line-num' data-line-number='313' href='#L313' id='L313'>
<i class="fa fa-link"></i>
313
</a>
<a class='diff-line-num' data-line-number='314' href='#L314' id='L314'>
<i class="fa fa-link"></i>
314
</a>
<a class='diff-line-num' data-line-number='315' href='#L315' id='L315'>
<i class="fa fa-link"></i>
315
</a>
<a class='diff-line-num' data-line-number='316' href='#L316' id='L316'>
<i class="fa fa-link"></i>
316
</a>
<a class='diff-line-num' data-line-number='317' href='#L317' id='L317'>
<i class="fa fa-link"></i>
317
</a>
<a class='diff-line-num' data-line-number='318' href='#L318' id='L318'>
<i class="fa fa-link"></i>
318
</a>
<a class='diff-line-num' data-line-number='319' href='#L319' id='L319'>
<i class="fa fa-link"></i>
319
</a>
<a class='diff-line-num' data-line-number='320' href='#L320' id='L320'>
<i class="fa fa-link"></i>
320
</a>
<a class='diff-line-num' data-line-number='321' href='#L321' id='L321'>
<i class="fa fa-link"></i>
321
</a>
<a class='diff-line-num' data-line-number='322' href='#L322' id='L322'>
<i class="fa fa-link"></i>
322
</a>
<a class='diff-line-num' data-line-number='323' href='#L323' id='L323'>
<i class="fa fa-link"></i>
323
</a>
<a class='diff-line-num' data-line-number='324' href='#L324' id='L324'>
<i class="fa fa-link"></i>
324
</a>
<a class='diff-line-num' data-line-number='325' href='#L325' id='L325'>
<i class="fa fa-link"></i>
325
</a>
<a class='diff-line-num' data-line-number='326' href='#L326' id='L326'>
<i class="fa fa-link"></i>
326
</a>
<a class='diff-line-num' data-line-number='327' href='#L327' id='L327'>
<i class="fa fa-link"></i>
327
</a>
<a class='diff-line-num' data-line-number='328' href='#L328' id='L328'>
<i class="fa fa-link"></i>
328
</a>
<a class='diff-line-num' data-line-number='329' href='#L329' id='L329'>
<i class="fa fa-link"></i>
329
</a>
<a class='diff-line-num' data-line-number='330' href='#L330' id='L330'>
<i class="fa fa-link"></i>
330
</a>
<a class='diff-line-num' data-line-number='331' href='#L331' id='L331'>
<i class="fa fa-link"></i>
331
</a>
<a class='diff-line-num' data-line-number='332' href='#L332' id='L332'>
<i class="fa fa-link"></i>
332
</a>
<a class='diff-line-num' data-line-number='333' href='#L333' id='L333'>
<i class="fa fa-link"></i>
333
</a>
<a class='diff-line-num' data-line-number='334' href='#L334' id='L334'>
<i class="fa fa-link"></i>
334
</a>
<a class='diff-line-num' data-line-number='335' href='#L335' id='L335'>
<i class="fa fa-link"></i>
335
</a>
<a class='diff-line-num' data-line-number='336' href='#L336' id='L336'>
<i class="fa fa-link"></i>
336
</a>
<a class='diff-line-num' data-line-number='337' href='#L337' id='L337'>
<i class="fa fa-link"></i>
337
</a>
<a class='diff-line-num' data-line-number='338' href='#L338' id='L338'>
<i class="fa fa-link"></i>
338
</a>
<a class='diff-line-num' data-line-number='339' href='#L339' id='L339'>
<i class="fa fa-link"></i>
339
</a>
<a class='diff-line-num' data-line-number='340' href='#L340' id='L340'>
<i class="fa fa-link"></i>
340
</a>
<a class='diff-line-num' data-line-number='341' href='#L341' id='L341'>
<i class="fa fa-link"></i>
341
</a>
<a class='diff-line-num' data-line-number='342' href='#L342' id='L342'>
<i class="fa fa-link"></i>
342
</a>
<a class='diff-line-num' data-line-number='343' href='#L343' id='L343'>
<i class="fa fa-link"></i>
343
</a>
<a class='diff-line-num' data-line-number='344' href='#L344' id='L344'>
<i class="fa fa-link"></i>
344
</a>
<a class='diff-line-num' data-line-number='345' href='#L345' id='L345'>
<i class="fa fa-link"></i>
345
</a>
<a class='diff-line-num' data-line-number='346' href='#L346' id='L346'>
<i class="fa fa-link"></i>
346
</a>
<a class='diff-line-num' data-line-number='347' href='#L347' id='L347'>
<i class="fa fa-link"></i>
347
</a>
<a class='diff-line-num' data-line-number='348' href='#L348' id='L348'>
<i class="fa fa-link"></i>
348
</a>
<a class='diff-line-num' data-line-number='349' href='#L349' id='L349'>
<i class="fa fa-link"></i>
349
</a>
<a class='diff-line-num' data-line-number='350' href='#L350' id='L350'>
<i class="fa fa-link"></i>
350
</a>
<a class='diff-line-num' data-line-number='351' href='#L351' id='L351'>
<i class="fa fa-link"></i>
351
</a>
<a class='diff-line-num' data-line-number='352' href='#L352' id='L352'>
<i class="fa fa-link"></i>
352
</a>
<a class='diff-line-num' data-line-number='353' href='#L353' id='L353'>
<i class="fa fa-link"></i>
353
</a>
<a class='diff-line-num' data-line-number='354' href='#L354' id='L354'>
<i class="fa fa-link"></i>
354
</a>
<a class='diff-line-num' data-line-number='355' href='#L355' id='L355'>
<i class="fa fa-link"></i>
355
</a>
<a class='diff-line-num' data-line-number='356' href='#L356' id='L356'>
<i class="fa fa-link"></i>
356
</a>
<a class='diff-line-num' data-line-number='357' href='#L357' id='L357'>
<i class="fa fa-link"></i>
357
</a>
<a class='diff-line-num' data-line-number='358' href='#L358' id='L358'>
<i class="fa fa-link"></i>
358
</a>
<a class='diff-line-num' data-line-number='359' href='#L359' id='L359'>
<i class="fa fa-link"></i>
359
</a>
<a class='diff-line-num' data-line-number='360' href='#L360' id='L360'>
<i class="fa fa-link"></i>
360
</a>
<a class='diff-line-num' data-line-number='361' href='#L361' id='L361'>
<i class="fa fa-link"></i>
361
</a>
<a class='diff-line-num' data-line-number='362' href='#L362' id='L362'>
<i class="fa fa-link"></i>
362
</a>
<a class='diff-line-num' data-line-number='363' href='#L363' id='L363'>
<i class="fa fa-link"></i>
363
</a>
<a class='diff-line-num' data-line-number='364' href='#L364' id='L364'>
<i class="fa fa-link"></i>
364
</a>
<a class='diff-line-num' data-line-number='365' href='#L365' id='L365'>
<i class="fa fa-link"></i>
365
</a>
<a class='diff-line-num' data-line-number='366' href='#L366' id='L366'>
<i class="fa fa-link"></i>
366
</a>
<a class='diff-line-num' data-line-number='367' href='#L367' id='L367'>
<i class="fa fa-link"></i>
367
</a>
<a class='diff-line-num' data-line-number='368' href='#L368' id='L368'>
<i class="fa fa-link"></i>
368
</a>
<a class='diff-line-num' data-line-number='369' href='#L369' id='L369'>
<i class="fa fa-link"></i>
369
</a>
<a class='diff-line-num' data-line-number='370' href='#L370' id='L370'>
<i class="fa fa-link"></i>
370
</a>
<a class='diff-line-num' data-line-number='371' href='#L371' id='L371'>
<i class="fa fa-link"></i>
371
</a>
<a class='diff-line-num' data-line-number='372' href='#L372' id='L372'>
<i class="fa fa-link"></i>
372
</a>
<a class='diff-line-num' data-line-number='373' href='#L373' id='L373'>
<i class="fa fa-link"></i>
373
</a>
<a class='diff-line-num' data-line-number='374' href='#L374' id='L374'>
<i class="fa fa-link"></i>
374
</a>
<a class='diff-line-num' data-line-number='375' href='#L375' id='L375'>
<i class="fa fa-link"></i>
375
</a>
<a class='diff-line-num' data-line-number='376' href='#L376' id='L376'>
<i class="fa fa-link"></i>
376
</a>
<a class='diff-line-num' data-line-number='377' href='#L377' id='L377'>
<i class="fa fa-link"></i>
377
</a>
<a class='diff-line-num' data-line-number='378' href='#L378' id='L378'>
<i class="fa fa-link"></i>
378
</a>
<a class='diff-line-num' data-line-number='379' href='#L379' id='L379'>
<i class="fa fa-link"></i>
379
</a>
<a class='diff-line-num' data-line-number='380' href='#L380' id='L380'>
<i class="fa fa-link"></i>
380
</a>
<a class='diff-line-num' data-line-number='381' href='#L381' id='L381'>
<i class="fa fa-link"></i>
381
</a>
<a class='diff-line-num' data-line-number='382' href='#L382' id='L382'>
<i class="fa fa-link"></i>
382
</a>
<a class='diff-line-num' data-line-number='383' href='#L383' id='L383'>
<i class="fa fa-link"></i>
383
</a>
<a class='diff-line-num' data-line-number='384' href='#L384' id='L384'>
<i class="fa fa-link"></i>
384
</a>
<a class='diff-line-num' data-line-number='385' href='#L385' id='L385'>
<i class="fa fa-link"></i>
385
</a>
<a class='diff-line-num' data-line-number='386' href='#L386' id='L386'>
<i class="fa fa-link"></i>
386
</a>
<a class='diff-line-num' data-line-number='387' href='#L387' id='L387'>
<i class="fa fa-link"></i>
387
</a>
<a class='diff-line-num' data-line-number='388' href='#L388' id='L388'>
<i class="fa fa-link"></i>
388
</a>
<a class='diff-line-num' data-line-number='389' href='#L389' id='L389'>
<i class="fa fa-link"></i>
389
</a>
<a class='diff-line-num' data-line-number='390' href='#L390' id='L390'>
<i class="fa fa-link"></i>
390
</a>
<a class='diff-line-num' data-line-number='391' href='#L391' id='L391'>
<i class="fa fa-link"></i>
391
</a>
<a class='diff-line-num' data-line-number='392' href='#L392' id='L392'>
<i class="fa fa-link"></i>
392
</a>
<a class='diff-line-num' data-line-number='393' href='#L393' id='L393'>
<i class="fa fa-link"></i>
393
</a>
<a class='diff-line-num' data-line-number='394' href='#L394' id='L394'>
<i class="fa fa-link"></i>
394
</a>
<a class='diff-line-num' data-line-number='395' href='#L395' id='L395'>
<i class="fa fa-link"></i>
395
</a>
<a class='diff-line-num' data-line-number='396' href='#L396' id='L396'>
<i class="fa fa-link"></i>
396
</a>
<a class='diff-line-num' data-line-number='397' href='#L397' id='L397'>
<i class="fa fa-link"></i>
397
</a>
<a class='diff-line-num' data-line-number='398' href='#L398' id='L398'>
<i class="fa fa-link"></i>
398
</a>
<a class='diff-line-num' data-line-number='399' href='#L399' id='L399'>
<i class="fa fa-link"></i>
399
</a>
<a class='diff-line-num' data-line-number='400' href='#L400' id='L400'>
<i class="fa fa-link"></i>
400
</a>
<a class='diff-line-num' data-line-number='401' href='#L401' id='L401'>
<i class="fa fa-link"></i>
401
</a>
<a class='diff-line-num' data-line-number='402' href='#L402' id='L402'>
<i class="fa fa-link"></i>
402
</a>
<a class='diff-line-num' data-line-number='403' href='#L403' id='L403'>
<i class="fa fa-link"></i>
403
</a>
<a class='diff-line-num' data-line-number='404' href='#L404' id='L404'>
<i class="fa fa-link"></i>
404
</a>
<a class='diff-line-num' data-line-number='405' href='#L405' id='L405'>
<i class="fa fa-link"></i>
405
</a>
<a class='diff-line-num' data-line-number='406' href='#L406' id='L406'>
<i class="fa fa-link"></i>
406
</a>
<a class='diff-line-num' data-line-number='407' href='#L407' id='L407'>
<i class="fa fa-link"></i>
407
</a>
<a class='diff-line-num' data-line-number='408' href='#L408' id='L408'>
<i class="fa fa-link"></i>
408
</a>
<a class='diff-line-num' data-line-number='409' href='#L409' id='L409'>
<i class="fa fa-link"></i>
409
</a>
<a class='diff-line-num' data-line-number='410' href='#L410' id='L410'>
<i class="fa fa-link"></i>
410
</a>
<a class='diff-line-num' data-line-number='411' href='#L411' id='L411'>
<i class="fa fa-link"></i>
411
</a>
<a class='diff-line-num' data-line-number='412' href='#L412' id='L412'>
<i class="fa fa-link"></i>
412
</a>
<a class='diff-line-num' data-line-number='413' href='#L413' id='L413'>
<i class="fa fa-link"></i>
413
</a>
<a class='diff-line-num' data-line-number='414' href='#L414' id='L414'>
<i class="fa fa-link"></i>
414
</a>
<a class='diff-line-num' data-line-number='415' href='#L415' id='L415'>
<i class="fa fa-link"></i>
415
</a>
<a class='diff-line-num' data-line-number='416' href='#L416' id='L416'>
<i class="fa fa-link"></i>
416
</a>
<a class='diff-line-num' data-line-number='417' href='#L417' id='L417'>
<i class="fa fa-link"></i>
417
</a>
<a class='diff-line-num' data-line-number='418' href='#L418' id='L418'>
<i class="fa fa-link"></i>
418
</a>
<a class='diff-line-num' data-line-number='419' href='#L419' id='L419'>
<i class="fa fa-link"></i>
419
</a>
<a class='diff-line-num' data-line-number='420' href='#L420' id='L420'>
<i class="fa fa-link"></i>
420
</a>
<a class='diff-line-num' data-line-number='421' href='#L421' id='L421'>
<i class="fa fa-link"></i>
421
</a>
<a class='diff-line-num' data-line-number='422' href='#L422' id='L422'>
<i class="fa fa-link"></i>
422
</a>
<a class='diff-line-num' data-line-number='423' href='#L423' id='L423'>
<i class="fa fa-link"></i>
423
</a>
<a class='diff-line-num' data-line-number='424' href='#L424' id='L424'>
<i class="fa fa-link"></i>
424
</a>
<a class='diff-line-num' data-line-number='425' href='#L425' id='L425'>
<i class="fa fa-link"></i>
425
</a>
<a class='diff-line-num' data-line-number='426' href='#L426' id='L426'>
<i class="fa fa-link"></i>
426
</a>
<a class='diff-line-num' data-line-number='427' href='#L427' id='L427'>
<i class="fa fa-link"></i>
427
</a>
<a class='diff-line-num' data-line-number='428' href='#L428' id='L428'>
<i class="fa fa-link"></i>
428
</a>
<a class='diff-line-num' data-line-number='429' href='#L429' id='L429'>
<i class="fa fa-link"></i>
429
</a>
<a class='diff-line-num' data-line-number='430' href='#L430' id='L430'>
<i class="fa fa-link"></i>
430
</a>
<a class='diff-line-num' data-line-number='431' href='#L431' id='L431'>
<i class="fa fa-link"></i>
431
</a>
<a class='diff-line-num' data-line-number='432' href='#L432' id='L432'>
<i class="fa fa-link"></i>
432
</a>
<a class='diff-line-num' data-line-number='433' href='#L433' id='L433'>
<i class="fa fa-link"></i>
433
</a>
<a class='diff-line-num' data-line-number='434' href='#L434' id='L434'>
<i class="fa fa-link"></i>
434
</a>
<a class='diff-line-num' data-line-number='435' href='#L435' id='L435'>
<i class="fa fa-link"></i>
435
</a>
<a class='diff-line-num' data-line-number='436' href='#L436' id='L436'>
<i class="fa fa-link"></i>
436
</a>
<a class='diff-line-num' data-line-number='437' href='#L437' id='L437'>
<i class="fa fa-link"></i>
437
</a>
<a class='diff-line-num' data-line-number='438' href='#L438' id='L438'>
<i class="fa fa-link"></i>
438
</a>
<a class='diff-line-num' data-line-number='439' href='#L439' id='L439'>
<i class="fa fa-link"></i>
439
</a>
<a class='diff-line-num' data-line-number='440' href='#L440' id='L440'>
<i class="fa fa-link"></i>
440
</a>
<a class='diff-line-num' data-line-number='441' href='#L441' id='L441'>
<i class="fa fa-link"></i>
441
</a>
<a class='diff-line-num' data-line-number='442' href='#L442' id='L442'>
<i class="fa fa-link"></i>
442
</a>
<a class='diff-line-num' data-line-number='443' href='#L443' id='L443'>
<i class="fa fa-link"></i>
443
</a>
<a class='diff-line-num' data-line-number='444' href='#L444' id='L444'>
<i class="fa fa-link"></i>
444
</a>
<a class='diff-line-num' data-line-number='445' href='#L445' id='L445'>
<i class="fa fa-link"></i>
445
</a>
<a class='diff-line-num' data-line-number='446' href='#L446' id='L446'>
<i class="fa fa-link"></i>
446
</a>
<a class='diff-line-num' data-line-number='447' href='#L447' id='L447'>
<i class="fa fa-link"></i>
447
</a>
<a class='diff-line-num' data-line-number='448' href='#L448' id='L448'>
<i class="fa fa-link"></i>
448
</a>
<a class='diff-line-num' data-line-number='449' href='#L449' id='L449'>
<i class="fa fa-link"></i>
449
</a>
<a class='diff-line-num' data-line-number='450' href='#L450' id='L450'>
<i class="fa fa-link"></i>
450
</a>
<a class='diff-line-num' data-line-number='451' href='#L451' id='L451'>
<i class="fa fa-link"></i>
451
</a>
<a class='diff-line-num' data-line-number='452' href='#L452' id='L452'>
<i class="fa fa-link"></i>
452
</a>
<a class='diff-line-num' data-line-number='453' href='#L453' id='L453'>
<i class="fa fa-link"></i>
453
</a>
<a class='diff-line-num' data-line-number='454' href='#L454' id='L454'>
<i class="fa fa-link"></i>
454
</a>
<a class='diff-line-num' data-line-number='455' href='#L455' id='L455'>
<i class="fa fa-link"></i>
455
</a>
<a class='diff-line-num' data-line-number='456' href='#L456' id='L456'>
<i class="fa fa-link"></i>
456
</a>
<a class='diff-line-num' data-line-number='457' href='#L457' id='L457'>
<i class="fa fa-link"></i>
457
</a>
<a class='diff-line-num' data-line-number='458' href='#L458' id='L458'>
<i class="fa fa-link"></i>
458
</a>
<a class='diff-line-num' data-line-number='459' href='#L459' id='L459'>
<i class="fa fa-link"></i>
459
</a>
<a class='diff-line-num' data-line-number='460' href='#L460' id='L460'>
<i class="fa fa-link"></i>
460
</a>
<a class='diff-line-num' data-line-number='461' href='#L461' id='L461'>
<i class="fa fa-link"></i>
461
</a>
<a class='diff-line-num' data-line-number='462' href='#L462' id='L462'>
<i class="fa fa-link"></i>
462
</a>
<a class='diff-line-num' data-line-number='463' href='#L463' id='L463'>
<i class="fa fa-link"></i>
463
</a>
<a class='diff-line-num' data-line-number='464' href='#L464' id='L464'>
<i class="fa fa-link"></i>
464
</a>
</div>
<div class='blob-content' data-blob-id='ec2c3ea37ebe4640e1e09aeeb2352a4fc130b7ce'>
<pre class="code highlight"><code><span id="LC1" class="line"><span class="n">function</span> <span class="n">ImaGIN_Electrode</span><span class="p">(</span><span class="n">S</span><span class="p">)</span></span>
<span id="LC2" class="line"><span class="o">%</span> <span class="n">Set</span> <span class="n">electrode</span> <span class="n">positions</span><span class="p">.</span></span>
<span id="LC3" class="line"></span>
<span id="LC4" class="line"><span class="o">%</span> <span class="o">-=============================================================================</span></span>
<span id="LC5" class="line"><span class="o">%</span> <span class="n">This</span> <span class="n">function</span> <span class="n">is</span> <span class="n">part</span> <span class="n">of</span> <span class="n">the</span> <span class="n">ImaGIN</span> <span class="n">software</span><span class="o">:</span> </span>
<span id="LC6" class="line"><span class="o">%</span> <span class="n">https</span><span class="o">:</span><span class="c1">//f-tract.eu/</span></span>
<span id="LC7" class="line"><span class="o">%</span></span>
<span id="LC8" class="line"><span class="o">%</span> <span class="n">This</span> <span class="n">software</span> <span class="n">is</span> <span class="n">distributed</span> <span class="n">under</span> <span class="n">the</span> <span class="n">terms</span> <span class="n">of</span> <span class="n">the</span> <span class="n">GNU</span> <span class="n">General</span> <span class="n">Public</span> <span class="n">License</span></span>
<span id="LC9" class="line"><span class="o">%</span> <span class="n">as</span> <span class="n">published</span> <span class="n">by</span> <span class="n">the</span> <span class="n">Free</span> <span class="n">Software</span> <span class="n">Foundation</span><span class="p">.</span> <span class="n">Further</span> <span class="n">details</span> <span class="n">on</span> <span class="n">the</span> <span class="n">GPLv3</span></span>
<span id="LC10" class="line"><span class="o">%</span> <span class="n">license</span> <span class="n">can</span> <span class="n">be</span> <span class="n">found</span> <span class="n">at</span> <span class="n">http</span><span class="o">:</span><span class="c1">//www.gnu.org/copyleft/gpl.html.</span></span>
<span id="LC11" class="line"><span class="o">%</span></span>
<span id="LC12" class="line"><span class="o">%</span> <span class="n">FOR</span> <span class="n">RESEARCH</span> <span class="n">PURPOSES</span> <span class="n">ONLY</span><span class="p">.</span> <span class="n">THE</span> <span class="n">SOFTWARE</span> <span class="n">IS</span> <span class="n">PROVIDED</span> <span class="s">&quot;AS IS,&quot;</span> <span class="n">AND</span> <span class="n">THE</span> <span class="n">AUTHORS</span></span>
<span id="LC13" class="line"><span class="o">%</span> <span class="n">DO</span> <span class="n">NOT</span> <span class="n">ASSUME</span> <span class="n">ANY</span> <span class="n">LIABILITY</span> <span class="n">OR</span> <span class="n">RESPONSIBILITY</span> <span class="n">FOR</span> <span class="n">ITS</span> <span class="n">USE</span> <span class="n">IN</span> <span class="n">ANY</span> <span class="n">CONTEXT</span><span class="p">.</span></span>
<span id="LC14" class="line"><span class="o">%</span></span>
<span id="LC15" class="line"><span class="o">%</span> <span class="n">Copyright</span> <span class="p">(</span><span class="n">c</span><span class="p">)</span> <span class="mi">2000</span><span class="o">-</span><span class="mi">2018</span> <span class="n">Inserm</span> <span class="n">U1216</span></span>
<span id="LC16" class="line"><span class="o">%</span> <span class="o">=============================================================================-</span></span>
<span id="LC17" class="line"><span class="o">%</span></span>
<span id="LC18" class="line"><span class="o">%</span> <span class="n">Authors</span><span class="o">:</span> <span class="n">Olivier</span> <span class="n">David</span><span class="p">,</span> <span class="n">Francois</span> <span class="n">Tadel</span></span>
<span id="LC19" class="line"></span>
<span id="LC20" class="line"><span class="o">%</span> <span class="o">%</span> <span class="n">Test</span> <span class="n">findChannel</span><span class="p">()</span></span>
<span id="LC21" class="line"><span class="o">%</span> <span class="n">ListCSV</span> <span class="o">=</span> <span class="p">{</span><span class="err">&#39;</span><span class="n">A01</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">A02</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">A03</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">A04</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">A05</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">A18</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Bp01</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Bp02</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">BP01</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">BP04</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Pp01</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">T101</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">T111</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">V101</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">V201</span><span class="err">&#39;</span><span class="p">};</span></span>
<span id="LC22" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">pp1</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">pp1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC23" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">A01</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">A01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC24" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">a01</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">a01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC25" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">a1</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>    <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">a1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC26" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">a</span> <span class="mi">1</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">a</span> <span class="mi">1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC27" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">A18</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">A18</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC28" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="sc">&#39;B&#39;</span><span class="err">&#39;</span><span class="mi">1</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="sc">&#39;B&#39;&#39;1&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC29" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="sc">&#39;b&#39;</span><span class="err">&#39;</span><span class="mo">01</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="sc">&#39;b&#39;</span><span class="err">&#39;</span><span class="mo">01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC30" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">Bp01</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Bp01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC31" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">Bp1</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Bp1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC32" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">bp1</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">bp1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC33" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span><span class="mi">1</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span><span class="mi">1</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC34" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span> <span class="mo">01</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span> <span class="mo">01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC35" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">BP01</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">BP01</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC36" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">Bp2</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Bp2</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC37" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="sc">&#39;B&#39;</span><span class="err">&#39;</span><span class="mi">2</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="sc">&#39;B&#39;&#39;2&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC38" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span> <span class="mi">2</span><span class="o">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">B</span><span class="p">,</span> <span class="mi">2</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC39" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">BP4</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">BP4</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC40" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">T111</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">T111</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC41" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">T11</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">T11</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC42" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">V101</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">V101</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC43" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">V11</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">V11</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC44" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">V101</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>  <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">V101</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC45" class="line"><span class="o">%</span> <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="nf">V21</span><span class="p">:</span> <span class="err">&#39;</span><span class="p">,</span>   <span class="n">ListCSV</span><span class="p">{</span><span class="n">findChannel</span><span class="p">(</span><span class="err">&#39;</span><span class="n">V21</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">ListCSV</span><span class="p">)}])</span></span>
<span id="LC46" class="line"><span class="o">%</span> <span class="k">return</span><span class="p">;</span></span>
<span id="LC47" class="line"></span>
<span id="LC48" class="line"><span class="o">%</span> <span class="n">Get</span> <span class="n">file</span> <span class="n">to</span> <span class="n">edit</span></span>
<span id="LC49" class="line"><span class="n">try</span></span>
<span id="LC50" class="line">    <span class="n">t</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">Fname</span><span class="p">;</span></span>
<span id="LC51" class="line"><span class="n">catch</span></span>
<span id="LC52" class="line">    <span class="n">t</span> <span class="o">=</span> <span class="n">spm_select</span><span class="p">(</span><span class="n">Inf</span><span class="p">,</span> <span class="err">&#39;\</span><span class="p">.</span><span class="n">mat</span><span class="err">$&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Select</span> <span class="n">data</span> <span class="n">file</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC53" class="line"><span class="n">end</span></span>
<span id="LC54" class="line"><span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">t</span><span class="p">)</span></span>
<span id="LC55" class="line">    <span class="k">return</span><span class="p">;</span></span>
<span id="LC56" class="line"><span class="n">end</span></span>
<span id="LC57" class="line"><span class="n">P</span> <span class="o">=</span> <span class="n">spm_str_manip</span><span class="p">(</span><span class="n">deblank</span><span class="p">(</span><span class="n">t</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span><span class="o">:</span><span class="p">)),</span><span class="sc">&#39;h&#39;</span><span class="p">);</span></span>
<span id="LC58" class="line"></span>
<span id="LC59" class="line"><span class="o">%</span> <span class="n">Read</span> <span class="n">sensor</span> <span class="n">names</span></span>
<span id="LC60" class="line"><span class="n">Position</span> <span class="o">=</span> <span class="p">[];</span></span>
<span id="LC61" class="line"><span class="n">try</span></span>
<span id="LC62" class="line">    <span class="n">Name</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">Name</span><span class="p">;</span></span>
<span id="LC63" class="line"><span class="n">catch</span></span>
<span id="LC64" class="line">    <span class="n">try</span></span>
<span id="LC65" class="line">        <span class="n">filename</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">filenameName</span><span class="p">;</span></span>
<span id="LC66" class="line">    <span class="n">catch</span></span>
<span id="LC67" class="line">        <span class="n">filename</span> <span class="o">=</span> <span class="n">spm_select</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="err">&#39;\</span><span class="p">.(</span><span class="n">csv</span><span class="o">|</span><span class="n">txt</span><span class="p">)</span><span class="err">$&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Select</span> <span class="n">txt</span> <span class="n">file</span> <span class="k">for</span> <span class="n">electrode</span> <span class="n">names</span><span class="err">&#39;</span><span class="p">,</span> <span class="p">{},</span> <span class="n">P</span><span class="p">);</span></span>
<span id="LC68" class="line">    <span class="n">end</span></span>
<span id="LC69" class="line">    <span class="o">%</span> <span class="n">Test</span> <span class="k">for</span> <span class="n">input</span> <span class="n">file</span></span>
<span id="LC70" class="line">    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">filename</span><span class="p">)</span> <span class="o">||</span> <span class="o">~</span><span class="n">exist</span><span class="p">(</span><span class="n">filename</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">file</span><span class="err">&#39;</span><span class="p">)</span></span>
<span id="LC71" class="line">        <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="n">ERROR</span><span class="o">:</span> <span class="n">Invalid</span> <span class="n">input</span> <span class="n">file</span> <span class="n">name</span><span class="p">.</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC72" class="line">        <span class="k">return</span></span>
<span id="LC73" class="line">    <span class="n">end</span></span>
<span id="LC74" class="line">    <span class="o">%</span> <span class="n">Detect</span> <span class="n">file</span> <span class="n">type</span><span class="o">:</span> <span class="n">_Name</span><span class="p">.</span><span class="n">txt</span><span class="o">/</span><span class="n">_Pos</span><span class="p">.</span><span class="n">txt</span> <span class="n">or</span> <span class="p">.</span><span class="n">csv</span></span>
<span id="LC75" class="line">    <span class="p">[</span><span class="n">fPath</span><span class="p">,</span> <span class="n">fBase</span><span class="p">,</span> <span class="nf">fExt</span><span class="p">]</span> <span class="o">=</span> <span class="n">fileparts</span><span class="p">(</span><span class="n">filename</span><span class="p">);</span></span>
<span id="LC76" class="line">    <span class="o">%</span> <span class="n">Read</span> <span class="n">file</span></span>
<span id="LC77" class="line">    <span class="k">switch</span> <span class="n">lower</span><span class="p">(</span><span class="n">fExt</span><span class="p">)</span></span>
<span id="LC78" class="line">        <span class="k">case</span> <span class="err">&#39;</span><span class="p">.</span><span class="n">txt</span><span class="err">&#39;</span></span>
<span id="LC79" class="line">            <span class="n">Name</span> <span class="o">=</span> <span class="n">readName</span><span class="p">(</span><span class="n">filename</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC80" class="line">        <span class="k">case</span> <span class="err">&#39;</span><span class="p">.</span><span class="n">csv</span><span class="err">&#39;</span></span>
<span id="LC81" class="line">            <span class="p">[</span><span class="n">Name</span><span class="p">,</span> <span class="nf">Position</span><span class="p">]</span> <span class="o">=</span> <span class="n">readCsv</span><span class="p">(</span><span class="n">filename</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC82" class="line">            <span class="k">if</span> <span class="n">all</span><span class="p">(</span><span class="n">isnan</span><span class="p">(</span><span class="n">Position</span><span class="p">))</span></span>
<span id="LC83" class="line">                <span class="n">error</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Anatomy</span> <span class="n">not</span> <span class="n">found</span><span class="p">:</span> <span class="n">possibly</span> <span class="n">patient</span> <span class="n">implanted</span> <span class="n">before</span> <span class="mi">2009</span><span class="p">.</span> <span class="n">Check</span> <span class="n">implantation</span> <span class="n">scheme</span><span class="o">!</span><span class="err">&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC84" class="line">            <span class="n">end</span></span>
<span id="LC85" class="line">        <span class="n">otherwise</span> </span>
<span id="LC86" class="line">            <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="n">ERROR</span><span class="p">:</span> <span class="n">Invalid</span> <span class="n">input</span> <span class="n">file</span> <span class="n">type</span><span class="p">.</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC87" class="line">            <span class="k">return</span></span>
<span id="LC88" class="line">    <span class="n">end</span></span>
<span id="LC89" class="line">    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">Name</span><span class="p">)</span></span>
<span id="LC90" class="line">        <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="n">ERROR</span><span class="o">:</span> <span class="n">No</span> <span class="n">contact</span> <span class="n">names</span> <span class="k">in</span> <span class="n">input</span> <span class="n">files</span><span class="p">.</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC91" class="line">        <span class="k">return</span></span>
<span id="LC92" class="line">    <span class="n">end</span></span>
<span id="LC93" class="line"><span class="n">end</span></span>
<span id="LC94" class="line"></span>
<span id="LC95" class="line"><span class="o">%</span> <span class="n">Read</span> <span class="n">sensor</span> <span class="n">positions</span></span>
<span id="LC96" class="line"><span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">Position</span><span class="p">)</span></span>
<span id="LC97" class="line">    <span class="n">try</span></span>
<span id="LC98" class="line">        <span class="n">Position</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">Position</span><span class="p">;</span></span>
<span id="LC99" class="line">    <span class="n">catch</span></span>
<span id="LC100" class="line">        <span class="n">try</span> </span>
<span id="LC101" class="line">            <span class="n">filename</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">filenamePos</span><span class="p">;</span></span>
<span id="LC102" class="line">            <span class="n">Position</span> <span class="o">=</span> <span class="n">load</span><span class="p">(</span><span class="n">filename</span><span class="p">);</span></span>
<span id="LC103" class="line">        <span class="n">catch</span></span>
<span id="LC104" class="line">            <span class="n">filename</span> <span class="o">=</span> <span class="n">spm_select</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="err">&#39;\</span><span class="p">.</span><span class="n">txt</span><span class="err">$&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Select</span> <span class="n">txt</span> <span class="n">file</span> <span class="k">for</span> <span class="n">electrode</span> <span class="n">positions</span><span class="err">&#39;</span><span class="p">,</span> <span class="p">{},</span> <span class="n">P</span><span class="p">);</span></span>
<span id="LC105" class="line">            <span class="n">Position</span> <span class="o">=</span> <span class="n">load</span><span class="p">(</span><span class="n">filename</span><span class="p">);</span></span>
<span id="LC106" class="line">        <span class="n">end</span></span>
<span id="LC107" class="line">    <span class="n">end</span></span>
<span id="LC108" class="line"><span class="n">end</span></span>
<span id="LC109" class="line"></span>
<span id="LC110" class="line"></span>
<span id="LC111" class="line"><span class="o">%</span> <span class="n">Set</span> <span class="n">positions</span></span>
<span id="LC112" class="line"><span class="n">chNotFound</span> <span class="o">=</span> <span class="p">{};</span></span>
<span id="LC113" class="line"><span class="n">chMatchLog</span> <span class="o">=</span> <span class="p">{};</span></span>
<span id="LC114" class="line"><span class="k">for</span> <span class="n">i0</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">size</span><span class="p">(</span><span class="n">t</span><span class="p">,</span><span class="mi">1</span><span class="p">)</span></span>
<span id="LC115" class="line">    <span class="n">T</span> <span class="o">=</span> <span class="n">deblank</span><span class="p">(</span><span class="n">t</span><span class="p">(</span><span class="n">i0</span><span class="p">,</span><span class="o">:</span><span class="p">));</span></span>
<span id="LC116" class="line">    <span class="o">%</span> <span class="n">Clone</span> <span class="n">file</span> <span class="k">if</span> <span class="n">requested</span> <span class="k">in</span> <span class="n">input</span></span>
<span id="LC117" class="line">    <span class="k">if</span> <span class="p">(</span><span class="n">nargin</span> <span class="o">&gt;=</span> <span class="mi">1</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="n">isfield</span><span class="p">(</span><span class="n">S</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">FileOut</span><span class="err">&#39;</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">S</span><span class="p">.</span><span class="n">FileOut</span><span class="p">)</span></span>
<span id="LC118" class="line">        <span class="n">D</span> <span class="o">=</span> <span class="n">spm_eeg_load</span><span class="p">(</span><span class="n">T</span><span class="p">);</span></span>
<span id="LC119" class="line">        <span class="n">D2</span> <span class="o">=</span> <span class="n">clone</span><span class="p">(</span><span class="n">D</span><span class="p">,</span> <span class="n">S</span><span class="p">.</span><span class="n">FileOut</span><span class="p">,</span> <span class="p">[</span><span class="n">D</span><span class="p">.</span><span class="n">nchannels</span> <span class="n">D</span><span class="p">.</span><span class="n">nsamples</span> <span class="n">D</span><span class="p">.</span><span class="nf">ntrials</span><span class="p">]);</span> </span>
<span id="LC120" class="line">        <span class="n">D2</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="o">:</span><span class="p">,</span><span class="o">:</span><span class="p">)</span> <span class="o">=</span> <span class="n">D</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="o">:</span><span class="p">,</span><span class="o">:</span><span class="p">);</span></span>
<span id="LC121" class="line">        <span class="n">save</span><span class="p">(</span><span class="n">D2</span><span class="p">);</span></span>
<span id="LC122" class="line">        <span class="n">SpmFile</span> <span class="o">=</span> <span class="n">S</span><span class="p">.</span><span class="n">FileOut</span><span class="p">;</span></span>
<span id="LC123" class="line">    <span class="k">else</span></span>
<span id="LC124" class="line">        <span class="n">SpmFile</span> <span class="o">=</span> <span class="n">T</span><span class="p">;</span></span>
<span id="LC125" class="line">    <span class="n">end</span></span>
<span id="LC126" class="line">    </span>
<span id="LC127" class="line">    <span class="o">%</span> <span class="n">Load</span> <span class="n">channel</span> <span class="n">names</span></span>
<span id="LC128" class="line">    <span class="n">SpmMat</span> <span class="o">=</span> <span class="n">load</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">);</span></span>
<span id="LC129" class="line">    <span class="n">Sensors</span> <span class="o">=</span> <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">sensors</span><span class="p">.</span><span class="n">eeg</span><span class="p">;</span></span>
<span id="LC130" class="line">    <span class="k">if</span> <span class="n">length</span><span class="p">(</span><span class="n">unique</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">))</span><span class="o">~=</span><span class="n">length</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">)</span></span>
<span id="LC131" class="line">        <span class="n">warning</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Repeated</span> <span class="n">label</span> <span class="k">in</span> <span class="n">the</span> <span class="n">file</span><span class="p">.</span><span class="err">&#39;</span><span class="p">)</span></span>
<span id="LC132" class="line">    <span class="n">end</span>    </span>
<span id="LC133" class="line">    <span class="o">%</span> <span class="n">Loop</span> <span class="n">on</span> <span class="n">all</span> <span class="n">channels</span> <span class="n">available</span> <span class="k">in</span> <span class="n">the</span> <span class="n">file</span></span>
<span id="LC134" class="line">    <span class="k">for</span> <span class="n">i1</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">length</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">)</span></span>
<span id="LC135" class="line">        <span class="n">sensLtmp</span> <span class="o">=</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">};</span></span>
<span id="LC136" class="line">        <span class="n">sensLtmp</span><span class="p">(</span><span class="n">ismember</span><span class="p">(</span><span class="kt">double</span><span class="p">(</span><span class="n">sensLtmp</span><span class="p">),[</span><span class="sc">&#39;,&#39;</span> <span class="sc">&#39;;&#39;</span> <span class="sc">&#39;-&#39;</span><span class="p">]))</span> <span class="o">=</span><span class="err">&#39;&#39;</span><span class="p">;</span></span>
<span id="LC137" class="line">        <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="o">=</span> <span class="n">sensLtmp</span><span class="p">;</span></span>
<span id="LC138" class="line">        <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">findChannel</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">},</span> <span class="n">Name</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">all_upper</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC139" class="line">        <span class="o">%</span> <span class="n">If</span> <span class="n">the</span> <span class="n">channel</span> <span class="n">was</span> <span class="n">already</span> <span class="n">found</span> <span class="k">in</span> <span class="n">the</span> <span class="n">list</span> <span class="n">before</span><span class="o">:</span> <span class="n">check</span> <span class="n">the</span> <span class="n">best</span> <span class="n">option</span> <span class="n">based</span> <span class="n">on</span> <span class="n">the</span> <span class="k">case</span></span>
<span id="LC140" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chMatchLog</span><span class="p">)</span></span>
<span id="LC141" class="line">            <span class="n">iPrevious</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmp</span><span class="p">(</span><span class="n">Name</span><span class="err">{</span><span class="n">iChanPos</span><span class="err">}</span><span class="p">,</span> <span class="n">chMatchLog</span><span class="p">(:,</span><span class="mi">2</span><span class="p">)));</span></span>
<span id="LC142" class="line">            <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iPrevious</span><span class="p">)</span> </span>
<span id="LC143" class="line">                <span class="o">%</span> <span class="n">If</span> <span class="n">the</span> <span class="n">new</span> <span class="n">channel</span> <span class="n">has</span> <span class="n">strictly</span> <span class="n">the</span> <span class="n">same</span> <span class="k">case</span><span class="p">,</span> <span class="n">or</span> <span class="k">if</span> <span class="n">it</span> <span class="n">corresponds</span> <span class="n">to</span> <span class="n">a</span> <span class="n">replaced</span> <span class="s">&quot;prime&quot;</span><span class="p">:</span> <span class="n">remove</span> <span class="n">the</span> <span class="n">previous</span> <span class="n">match</span></span>
<span id="LC144" class="line">                <span class="k">if</span> <span class="n">isequal</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">},</span> <span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">})</span> <span class="o">||</span> <span class="p">...</span></span>
<span id="LC145" class="line">                   <span class="p">(</span><span class="n">any</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="o">==</span> <span class="err">&#39;&#39;&#39;&#39;</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="n">strcmp</span><span class="p">(</span><span class="n">strrep</span><span class="p">(</span><span class="n">upper</span><span class="p">(</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}),</span> <span class="err">&#39;&#39;&#39;&#39;</span><span class="p">,</span> <span class="sc">&#39;p&#39;</span><span class="p">),</span> <span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">}))</span></span>
<span id="LC146" class="line">                    <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="nf">WARNING</span><span class="p">:</span> <span class="n">Channel</span> <span class="n">name</span> <span class="nf">conflict</span><span class="p">:</span> <span class="err">&#39;</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="err">&#39;</span> <span class="n">matched</span> <span class="n">with</span> <span class="err">&#39;</span> <span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">}</span> <span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">iPrevious</span><span class="p">,</span><span class="mi">1</span><span class="p">}</span> <span class="err">&#39;</span> <span class="n">discarded</span><span class="p">.</span><span class="err">&#39;</span><span class="p">]);</span></span>
<span id="LC147" class="line">                    <span class="n">chMatchLog</span><span class="p">(</span><span class="n">iPrevious</span><span class="p">,</span><span class="o">:</span><span class="p">)</span> <span class="o">=</span> <span class="p">[];</span>                </span>
<span id="LC148" class="line">                <span class="o">%</span> <span class="n">If</span> <span class="n">both</span> <span class="n">labels</span> <span class="n">are</span> <span class="n">the</span> <span class="n">same</span> <span class="n">regardless</span> <span class="n">of</span> <span class="n">the</span> <span class="k">case</span><span class="p">,</span> <span class="n">raises</span> <span class="n">an</span> <span class="n">error</span><span class="p">.</span></span>
<span id="LC149" class="line">                <span class="n">elseif</span> <span class="n">strcmpi</span><span class="p">(</span><span class="n">chMatchLog</span><span class="err">{</span><span class="n">iPrevious</span><span class="p">,</span><span class="mi">1</span><span class="err">}</span><span class="p">,</span><span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="err">{</span><span class="n">i1</span><span class="err">}</span><span class="p">)</span></span>
<span id="LC150" class="line">                    <span class="n">error</span><span class="p">([</span><span class="err">&#39;</span><span class="n">Duplicated</span> <span class="n">label</span> <span class="nf">found</span><span class="p">:</span> <span class="err">&#39;</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}])</span><span class="err">;</span>                                        </span>
<span id="LC151" class="line">                <span class="k">else</span></span>
<span id="LC152" class="line">                    <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="nf">WARNING</span><span class="p">:</span> <span class="n">Channel</span> <span class="n">name</span> <span class="nf">conflict</span><span class="p">:</span> <span class="err">&#39;</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">iPrevious</span><span class="p">,</span><span class="mi">1</span><span class="p">}</span> <span class="err">&#39;</span> <span class="n">matched</span> <span class="n">with</span> <span class="err">&#39;</span> <span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">}</span> <span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="err">&#39;</span> <span class="n">discarded</span><span class="p">.</span><span class="err">&#39;</span><span class="p">])</span><span class="err">;</span></span>
<span id="LC153" class="line">                    <span class="n">iChanPos</span> <span class="o">=</span> <span class="p">[]</span><span class="err">;</span></span>
<span id="LC154" class="line">                <span class="n">end</span></span>
<span id="LC155" class="line">            <span class="n">end</span></span>
<span id="LC156" class="line">        <span class="n">end</span></span>
<span id="LC157" class="line">        <span class="o">%</span> <span class="n">If</span> <span class="n">channel</span> <span class="n">was</span> <span class="n">found</span></span>
<span id="LC158" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC159" class="line">            <span class="o">%</span> <span class="n">Copy</span> <span class="n">position</span></span>
<span id="LC160" class="line">            <span class="n">Sensors</span><span class="p">.</span><span class="n">elecpos</span><span class="p">(</span><span class="n">i1</span><span class="p">,:)</span> <span class="o">=</span> <span class="n">Position</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">,</span><span class="o">:</span><span class="p">);</span></span>
<span id="LC161" class="line">            <span class="n">Sensors</span><span class="p">.</span><span class="n">chanpos</span><span class="p">(</span><span class="n">i1</span><span class="p">,</span><span class="o">:</span><span class="p">)</span> <span class="o">=</span> <span class="n">Position</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">,</span><span class="o">:</span><span class="p">);</span></span>
<span id="LC162" class="line">            <span class="o">%</span> <span class="n">Copy</span> <span class="n">channel</span> <span class="n">name</span> <span class="n">from</span> <span class="n">input</span> <span class="n">name</span> <span class="n">file</span> <span class="p">(</span><span class="n">ADDED</span> <span class="n">BY</span> <span class="n">FT</span> <span class="mi">5</span><span class="o">-</span><span class="n">Oct</span><span class="o">-</span><span class="mi">2018</span><span class="p">)</span></span>
<span id="LC163" class="line">            <span class="n">chMatchLog</span><span class="p">{</span><span class="n">end</span><span class="o">+</span><span class="mi">1</span><span class="p">,</span><span class="mi">1</span><span class="p">}</span> <span class="o">=</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">};</span></span>
<span id="LC164" class="line">            <span class="n">chMatchLog</span><span class="p">{</span><span class="n">end</span><span class="p">,</span><span class="mi">2</span><span class="p">}</span> <span class="o">=</span> <span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">};</span></span>
<span id="LC165" class="line">            <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">Name</span><span class="p">{</span><span class="n">iChanPos</span><span class="p">},</span><span class="sc">&#39;-&#39;</span><span class="p">,</span><span class="err">&#39;&#39;</span><span class="p">);</span></span>
<span id="LC166" class="line">        <span class="k">else</span></span>
<span id="LC167" class="line">            <span class="n">disp</span><span class="p">([</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="nf">WARNING</span><span class="p">:</span> <span class="err">&#39;</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="err">&#39;</span> <span class="n">not</span> <span class="n">assigned</span><span class="err">&#39;</span><span class="p">]);</span></span>
<span id="LC168" class="line">            <span class="n">chNotFound</span><span class="p">{</span><span class="n">end</span><span class="o">+</span><span class="mi">1</span><span class="p">}</span> <span class="o">=</span> <span class="n">Sensors</span><span class="p">.</span><span class="n">label</span><span class="p">{</span><span class="n">i1</span><span class="p">};</span></span>
<span id="LC169" class="line">        <span class="n">end</span></span>
<span id="LC170" class="line">    <span class="n">end</span></span>
<span id="LC171" class="line">    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">chMatchLog</span><span class="p">)</span></span>
<span id="LC172" class="line">        <span class="n">chMatchLog</span> <span class="o">=</span> <span class="n">repmat</span><span class="p">({</span><span class="err">&#39;&#39;</span><span class="p">},</span><span class="mi">1</span><span class="p">,</span><span class="mi">2</span><span class="p">);</span>        </span>
<span id="LC173" class="line">    <span class="n">end</span></span>
<span id="LC174" class="line">    <span class="o">%</span> <span class="n">Electrodes</span> <span class="n">present</span> <span class="k">in</span> <span class="n">the</span> <span class="n">CSV</span> <span class="n">but</span> <span class="n">not</span> <span class="n">found</span> <span class="k">in</span> <span class="n">the</span> <span class="n">SEEG</span> <span class="n">recordings</span></span>
<span id="LC175" class="line">    <span class="p">[</span><span class="o">~</span><span class="p">,</span> <span class="o">~</span><span class="p">,</span> <span class="nf">chTagsCSV</span><span class="p">]</span> <span class="o">=</span> <span class="n">ImaGIN_select_channels</span><span class="p">(</span><span class="n">Name</span><span class="p">);</span></span>
<span id="LC176" class="line">    <span class="p">[</span><span class="o">~</span><span class="p">,</span> <span class="o">~</span><span class="p">,</span> <span class="nf">chTagsSEEG</span><span class="p">]</span> <span class="o">=</span> <span class="n">ImaGIN_select_channels</span><span class="p">(</span><span class="n">chMatchLog</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="mi">2</span><span class="p">)</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC177" class="line">    <span class="n">elecUnused</span> <span class="o">=</span> <span class="n">setdiff</span><span class="p">(</span><span class="n">unique</span><span class="p">(</span><span class="n">chTagsCSV</span><span class="p">),</span> <span class="n">unique</span><span class="p">(</span><span class="n">chTagsSEEG</span><span class="p">));</span></span>
<span id="LC178" class="line">    </span>
<span id="LC179" class="line">    <span class="o">%</span> <span class="n">Add</span> <span class="n">entries</span> <span class="k">in</span> <span class="n">a</span> <span class="n">NEW</span> <span class="n">log</span> <span class="n">file</span></span>
<span id="LC180" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chMatchLog</span><span class="p">)</span> </span>
<span id="LC181" class="line">        <span class="n">ImaGIN_save_log</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Positions</span> <span class="n">added</span> <span class="k">for</span> <span class="n">channels</span><span class="o">:</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chMatchLog</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="mi">1</span><span class="p">));</span></span>
<span id="LC182" class="line">    <span class="n">end</span></span>
<span id="LC183" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chNotFound</span><span class="p">)</span> </span>
<span id="LC184" class="line">        <span class="n">ImaGIN_save_log</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Unmatched</span> <span class="n">SEEG</span> <span class="n">channels</span><span class="o">:</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chNotFound</span><span class="p">);</span></span>
<span id="LC185" class="line">        <span class="o">%</span> <span class="n">save</span> <span class="n">unmatched</span> <span class="n">channels</span> <span class="k">in</span> <span class="n">separate</span> <span class="p">.</span><span class="n">txt</span> <span class="n">file</span></span>
<span id="LC186" class="line">        <span class="p">[</span><span class="n">unPath</span><span class="p">,</span> <span class="n">unFile</span><span class="p">,</span> <span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">fileparts</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">);</span></span>
<span id="LC187" class="line">        <span class="k">if</span> <span class="n">numel</span><span class="p">(</span><span class="n">chNotFound</span><span class="p">)</span> <span class="o">&gt;</span> <span class="n">numel</span><span class="p">(</span><span class="n">cell2mat</span><span class="p">(</span><span class="n">regexpi</span><span class="p">(</span><span class="n">chNotFound</span><span class="p">,</span><span class="err">&#39;</span><span class="n">ecg</span><span class="err">&#39;</span><span class="p">)))</span></span>
<span id="LC188" class="line">            <span class="n">try</span></span>
<span id="LC189" class="line">                <span class="n">unmatchedFileName</span>  <span class="o">=</span> <span class="n">fullfile</span><span class="p">(</span><span class="n">unPath</span><span class="p">,</span> <span class="p">[</span><span class="err">&#39;</span><span class="n">contactsunmatched_</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">unFile</span><span class="p">,</span><span class="err">&#39;</span><span class="p">.</span><span class="n">txt</span><span class="err">&#39;</span><span class="p">]);</span></span>
<span id="LC190" class="line">                <span class="n">fid</span> <span class="o">=</span> <span class="n">fopen</span><span class="p">(</span><span class="n">unmatchedFileName</span><span class="p">,</span><span class="sc">&#39;w&#39;</span><span class="p">);</span></span>
<span id="LC191" class="line">                <span class="k">for</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">length</span><span class="p">(</span><span class="n">chNotFound</span><span class="p">)</span></span>
<span id="LC192" class="line">                    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">regexpi</span><span class="p">(</span><span class="n">chNotFound</span><span class="p">{</span><span class="n">i</span><span class="p">},</span><span class="err">&#39;</span><span class="n">ecg</span><span class="err">&#39;</span><span class="p">))</span></span>
<span id="LC193" class="line">                        <span class="n">fprintf</span><span class="p">(</span><span class="n">fid</span><span class="p">,</span><span class="err">&#39;</span><span class="o">%</span><span class="n">s</span><span class="err">\</span><span class="n">n</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chNotFound</span><span class="p">{</span><span class="n">i</span><span class="p">});</span></span>
<span id="LC194" class="line">                    <span class="n">end</span></span>
<span id="LC195" class="line">                <span class="n">end</span></span>
<span id="LC196" class="line">                <span class="n">fclose</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC197" class="line">            <span class="n">catch</span></span>
<span id="LC198" class="line">                <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Unmatched</span> <span class="n">channels</span> <span class="n">names</span> <span class="n">SEEG</span><span class="o">-</span><span class="n">CSV</span> <span class="n">not</span> <span class="n">saved</span><span class="p">.</span><span class="err">&#39;</span><span class="p">)</span></span>
<span id="LC199" class="line">            <span class="n">end</span></span>
<span id="LC200" class="line">        <span class="n">end</span></span>
<span id="LC201" class="line">    <span class="n">end</span></span>
<span id="LC202" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">elecUnused</span><span class="p">)</span> </span>
<span id="LC203" class="line">        <span class="n">ImaGIN_save_log</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">Unmatched</span> <span class="n">CSV</span> <span class="n">electrodes</span><span class="o">:</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">elecUnused</span><span class="p">);</span></span>
<span id="LC204" class="line">    <span class="n">end</span></span>
<span id="LC205" class="line">    </span>
<span id="LC206" class="line">    <span class="o">%</span> <span class="n">Match</span> <span class="n">file</span><span class="o">:</span> <span class="n">Correspondance</span> <span class="n">between</span> <span class="n">SEEG</span><span class="o">-</span><span class="n">CSV</span><span class="o">-</span><span class="n">LENA</span> <span class="n">conventions</span></span>
<span id="LC207" class="line">    <span class="k">if</span> <span class="n">isfield</span><span class="p">(</span><span class="n">S</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">FileTxtOut</span><span class="err">&#39;</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">S</span><span class="p">.</span><span class="n">FileTxtOut</span><span class="p">)</span></span>
<span id="LC208" class="line">        <span class="n">try</span></span>
<span id="LC209" class="line">            <span class="n">fid</span> <span class="o">=</span> <span class="n">fopen</span><span class="p">(</span><span class="n">S</span><span class="p">.</span><span class="n">FileTxtOut</span><span class="p">,</span><span class="sc">&#39;w&#39;</span><span class="p">);</span></span>
<span id="LC210" class="line">            <span class="n">fprintf</span><span class="p">(</span><span class="n">fid</span><span class="p">,</span><span class="err">&#39;</span><span class="n">SEEG</span><span class="p">,</span><span class="n">CSV</span><span class="err">\</span><span class="n">n</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC211" class="line">            <span class="k">for</span> <span class="n">i</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">size</span><span class="p">(</span><span class="n">chMatchLog</span><span class="p">,</span><span class="mi">1</span><span class="p">)</span></span>
<span id="LC212" class="line">                <span class="n">fprintf</span><span class="p">(</span><span class="n">fid</span><span class="p">,</span><span class="err">&#39;</span><span class="o">%</span><span class="n">s</span><span class="p">,</span><span class="o">%</span><span class="n">s</span><span class="err">\</span><span class="n">n</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">i</span><span class="p">,</span><span class="mi">1</span><span class="p">},</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">i</span><span class="p">,</span><span class="mi">2</span><span class="p">});</span></span>
<span id="LC213" class="line">            <span class="n">end</span></span>
<span id="LC214" class="line">            <span class="n">fclose</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC215" class="line">        <span class="n">catch</span></span>
<span id="LC216" class="line">            <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Log</span> <span class="n">with</span> <span class="n">matched</span> <span class="n">channels</span> <span class="n">names</span> <span class="n">SEEG</span><span class="o">-</span><span class="n">CSV</span> <span class="n">not</span> <span class="n">saved</span><span class="p">.</span><span class="err">&#39;</span><span class="p">)</span></span>
<span id="LC217" class="line">        <span class="n">end</span></span>
<span id="LC218" class="line">    <span class="n">end</span></span>
<span id="LC219" class="line">    </span>
<span id="LC220" class="line">    <span class="o">%</span> <span class="n">Replace</span> <span class="n">channel</span> <span class="n">definition</span> <span class="k">in</span> <span class="n">input</span> <span class="p">.</span><span class="n">mat</span><span class="o">/</span><span class="p">.</span><span class="n">dat</span> <span class="n">file</span></span>
<span id="LC221" class="line">    <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">sensors</span><span class="p">.</span><span class="n">eeg</span> <span class="o">=</span> <span class="n">Sensors</span><span class="p">;</span></span>
<span id="LC222" class="line">    </span>
<span id="LC223" class="line">    <span class="o">%</span> <span class="n">Replace</span> <span class="n">labels</span> <span class="k">in</span> <span class="n">D</span><span class="p">.</span><span class="n">channels</span></span>
<span id="LC224" class="line">    <span class="k">for</span> <span class="n">iChan</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">length</span><span class="p">(</span><span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">channels</span><span class="p">)</span></span>
<span id="LC225" class="line">        <span class="n">spmLtmp</span> <span class="o">=</span> <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">channels</span><span class="p">(</span><span class="n">iChan</span><span class="p">).</span><span class="n">label</span><span class="p">;</span></span>
<span id="LC226" class="line">        <span class="n">spmLtmp</span><span class="p">(</span><span class="n">ismember</span><span class="p">(</span><span class="kt">double</span><span class="p">(</span><span class="n">spmLtmp</span><span class="p">),[</span><span class="sc">&#39;,&#39;</span> <span class="sc">&#39;;&#39;</span> <span class="sc">&#39;-&#39;</span><span class="p">]))</span> <span class="o">=</span><span class="err">&#39;&#39;</span><span class="p">;</span></span>
<span id="LC227" class="line">        <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">channels</span><span class="p">(</span><span class="n">iChan</span><span class="p">).</span><span class="n">label</span> <span class="o">=</span> <span class="n">spmLtmp</span><span class="p">;</span></span>
<span id="LC228" class="line">        </span>
<span id="LC229" class="line">        <span class="n">iChanMatch</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">channels</span><span class="p">(</span><span class="n">iChan</span><span class="p">).</span><span class="n">label</span><span class="p">,</span> <span class="n">chMatchLog</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="mi">1</span><span class="p">)));</span></span>
<span id="LC230" class="line">        <span class="k">if</span> <span class="p">(</span><span class="n">length</span><span class="p">(</span><span class="n">iChanMatch</span><span class="p">)</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span></span>
<span id="LC231" class="line">            <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">channels</span><span class="p">(</span><span class="n">iChan</span><span class="p">).</span><span class="n">label</span> <span class="o">=</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">iChanMatch</span><span class="p">,</span><span class="mi">2</span><span class="p">};</span></span>
<span id="LC232" class="line">        <span class="n">end</span></span>
<span id="LC233" class="line">    <span class="n">end</span></span>
<span id="LC234" class="line">    </span>
<span id="LC235" class="line">    <span class="o">%</span> <span class="n">Replace</span> <span class="n">labels</span> <span class="k">in</span> <span class="n">events</span></span>
<span id="LC236" class="line">    <span class="k">for</span> <span class="n">iEvt</span> <span class="o">=</span> <span class="mi">1</span><span class="o">:</span><span class="n">length</span><span class="p">(</span><span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">trials</span><span class="p">.</span><span class="n">events</span><span class="p">)</span></span>
<span id="LC237" class="line">        <span class="o">%</span> <span class="n">Get</span> <span class="n">channel</span> <span class="n">name</span> <span class="n">from</span> <span class="n">event</span> <span class="n">name</span></span>
<span id="LC238" class="line">        <span class="n">EvtName</span> <span class="o">=</span> <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">trials</span><span class="p">.</span><span class="n">events</span><span class="p">(</span><span class="n">iEvt</span><span class="p">).</span><span class="n">type</span><span class="p">;</span></span>
<span id="LC239" class="line">        <span class="p">[</span><span class="n">chLabel1</span><span class="p">,</span> <span class="n">chLabel2</span><span class="p">,</span> <span class="n">noteNameNew</span><span class="p">,</span> <span class="n">chInd1</span><span class="p">,</span> <span class="nf">chInd2</span><span class="p">]</span> <span class="o">=</span> <span class="n">ImaGIN_CleanEventName</span><span class="p">(</span><span class="n">EvtName</span><span class="p">);</span>  </span>
<span id="LC240" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC241" class="line">            <span class="k">if</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">&lt;</span>  <span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC242" class="line">                <span class="k">if</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">+</span><span class="mi">1</span> <span class="o">~=</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC243" class="line">                    <span class="n">chInd2tmp</span> <span class="o">=</span> <span class="n">num2str</span><span class="p">(</span><span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">+</span> <span class="mi">1</span><span class="p">);</span></span>
<span id="LC244" class="line">                    <span class="k">if</span> <span class="n">numel</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">~=</span> <span class="n">numel</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC245" class="line">                        <span class="n">chInd2tmp</span> <span class="o">=</span> <span class="n">strcat</span><span class="p">(</span><span class="sc">&#39;0&#39;</span><span class="p">,</span><span class="n">chInd2tmp</span><span class="p">);</span></span>
<span id="LC246" class="line">                    <span class="n">end</span></span>
<span id="LC247" class="line">                    <span class="n">chLabel2</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel2</span><span class="p">,</span><span class="n">chInd2</span><span class="p">,</span><span class="n">chInd2tmp</span><span class="p">);</span></span>
<span id="LC248" class="line">                    <span class="n">chInd2</span>   <span class="o">=</span> <span class="n">chInd2tmp</span><span class="p">;</span></span>
<span id="LC249" class="line">                <span class="n">end</span></span>
<span id="LC250" class="line">            <span class="n">elseif</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">&gt;</span>  <span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC251" class="line">                <span class="k">if</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span> <span class="o">~=</span> <span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">1</span></span>
<span id="LC252" class="line">                    <span class="n">chInd1tmp</span> <span class="o">=</span> <span class="n">num2str</span><span class="p">(</span><span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span> <span class="o">+</span> <span class="mi">1</span><span class="p">);</span></span>
<span id="LC253" class="line">                    <span class="k">if</span> <span class="n">numel</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span> <span class="o">~=</span> <span class="n">numel</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span></span>
<span id="LC254" class="line">                        <span class="n">chInd1tmp</span> <span class="o">=</span> <span class="n">strcat</span><span class="p">(</span><span class="sc">&#39;0&#39;</span><span class="p">,</span><span class="n">chInd1</span><span class="p">);</span></span>
<span id="LC255" class="line">                    <span class="n">end</span></span>
<span id="LC256" class="line">                    <span class="n">chLabel1</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel1</span><span class="p">,</span><span class="n">chInd1</span><span class="p">,</span><span class="n">chInd1tmp</span><span class="p">);</span></span>
<span id="LC257" class="line">                    <span class="n">chInd1</span>   <span class="o">=</span> <span class="n">chInd1tmp</span><span class="p">;</span></span>
<span id="LC258" class="line">                <span class="n">end</span> </span>
<span id="LC259" class="line">            <span class="n">end</span></span>
<span id="LC260" class="line">        <span class="n">end</span></span>
<span id="LC261" class="line">        <span class="o">%</span> <span class="n">Replace</span> <span class="n">with</span> <span class="n">matching</span> <span class="n">CSV</span> <span class="n">name</span></span>
<span id="LC262" class="line">        <span class="o">%</span> <span class="n">This</span> <span class="n">section</span> <span class="n">now</span> <span class="n">uses</span> <span class="s">&quot;Name&quot;</span> <span class="n">extracted</span> <span class="n">from</span> <span class="n">the</span> <span class="n">csv</span> <span class="n">to</span> <span class="n">consider</span> <span class="n">all</span> <span class="n">electrode</span> <span class="n">labels</span> <span class="n">instead</span> <span class="n">of</span> <span class="n">the</span> <span class="n">ones</span> <span class="k">in</span> <span class="s">&quot;chMatchLog&quot;</span> <span class="n">as</span> <span class="n">it</span> <span class="n">used</span> <span class="n">to</span> <span class="k">do</span> <span class="n">because</span> <span class="n">it</span> <span class="n">only</span> <span class="n">considered</span> <span class="n">electrodes</span> <span class="n">which</span> <span class="n">also</span> <span class="n">recorded</span><span class="p">.</span></span>
<span id="LC263" class="line">        <span class="o">%</span> <span class="n">Note</span> <span class="n">by</span> <span class="n">A</span><span class="p">.</span> <span class="n">Boyer</span> <span class="o">-</span> <span class="mi">13</span><span class="o">/</span><span class="mo">02</span><span class="o">/</span><span class="mi">2020</span></span>
<span id="LC264" class="line">        <span class="n">elec_labels</span> <span class="o">=</span> <span class="n">Name</span><span class="err">&#39;</span><span class="p">;</span></span>
<span id="LC265" class="line">        <span class="n">elec_labels_no_primes</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">elec_labels</span><span class="p">,</span><span class="sc">&#39;p&#39;</span><span class="p">,</span><span class="err">&#39;&#39;&#39;&#39;</span><span class="p">);</span> </span>
<span id="LC266" class="line">        <span class="n">csv_all_electrodes</span> <span class="o">=</span> <span class="p">[</span><span class="n">elec_labels</span> <span class="nf">elec_labels_no_primes</span><span class="p">];</span>      </span>
<span id="LC267" class="line">        </span>
<span id="LC268" class="line">        <span class="p">[</span><span class="n">iChanMatch1</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">chLabel1</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span></span>
<span id="LC269" class="line">        <span class="p">[</span><span class="n">iChanMatch2</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">chLabel2</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span>        </span>
<span id="LC270" class="line">        <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)</span></span>
<span id="LC271" class="line">            <span class="n">tmpchLabel1</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel1</span><span class="p">,</span> <span class="n">chInd1</span><span class="p">,</span> <span class="n">num2str</span><span class="p">(</span><span class="n">str2double</span><span class="p">(</span><span class="n">chInd1</span><span class="p">)));</span></span>
<span id="LC272" class="line">            <span class="p">[</span><span class="n">iChanMatch1</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">tmpchLabel1</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span></span>
<span id="LC273" class="line">            <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span></span>
<span id="LC274" class="line">                <span class="n">tmpchLabel1</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel1</span><span class="p">,</span> <span class="n">chInd1</span><span class="p">,</span> <span class="p">[</span><span class="sc">&#39;0&#39;</span> <span class="nf">chInd1</span><span class="p">]);</span></span>
<span id="LC275" class="line">                <span class="p">[</span><span class="n">iChanMatch1</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">tmpchLabel1</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span></span>
<span id="LC276" class="line">            <span class="n">end</span>   </span>
<span id="LC277" class="line">        <span class="n">end</span>        </span>
<span id="LC278" class="line">        <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)</span></span>
<span id="LC279" class="line">            <span class="n">tmpchLabel2</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel2</span><span class="p">,</span> <span class="n">chInd2</span><span class="p">,</span> <span class="n">num2str</span><span class="p">(</span><span class="n">str2double</span><span class="p">(</span><span class="n">chInd2</span><span class="p">)));</span></span>
<span id="LC280" class="line">            <span class="p">[</span><span class="n">iChanMatch2</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">tmpchLabel2</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span></span>
<span id="LC281" class="line">            <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span></span>
<span id="LC282" class="line">                <span class="n">tmpchLabel2</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">chLabel2</span><span class="p">,</span> <span class="n">chInd2</span><span class="p">,</span> <span class="p">[</span><span class="sc">&#39;0&#39;</span> <span class="nf">chInd2</span><span class="p">]);</span></span>
<span id="LC283" class="line">                <span class="p">[</span><span class="n">iChanMatch2</span><span class="p">,</span><span class="o">~</span><span class="p">]</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">tmpchLabel2</span><span class="p">,</span> <span class="n">csv_all_electrodes</span><span class="p">));</span></span>
<span id="LC284" class="line">            <span class="n">end</span></span>
<span id="LC285" class="line">        <span class="n">end</span></span>
<span id="LC286" class="line">        </span>
<span id="LC287" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span></span>
<span id="LC288" class="line">            <span class="k">if</span> <span class="n">sum</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="o">==</span><span class="n">iChanMatch1</span><span class="p">(</span><span class="mi">1</span><span class="p">))</span> <span class="o">==</span> <span class="n">numel</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span></span>
<span id="LC289" class="line">                <span class="n">iChanMatch1</span> <span class="o">=</span> <span class="n">iChanMatch1</span><span class="p">(</span><span class="mi">1</span><span class="p">);</span></span>
<span id="LC290" class="line">            <span class="k">else</span></span>
<span id="LC291" class="line">                <span class="n">warning</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Same</span> <span class="n">label</span> <span class="n">found</span> <span class="k">for</span> <span class="n">multiple</span> <span class="n">electrodes</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC292" class="line">            <span class="n">end</span></span>
<span id="LC293" class="line">        <span class="n">end</span></span>
<span id="LC294" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span></span>
<span id="LC295" class="line">            <span class="k">if</span> <span class="n">sum</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="o">==</span><span class="n">iChanMatch2</span><span class="p">(</span><span class="mi">1</span><span class="p">))</span> <span class="o">==</span> <span class="n">numel</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span></span>
<span id="LC296" class="line">                <span class="n">iChanMatch2</span> <span class="o">=</span> <span class="n">iChanMatch2</span><span class="p">(</span><span class="mi">1</span><span class="p">);</span></span>
<span id="LC297" class="line">            <span class="k">else</span></span>
<span id="LC298" class="line">                <span class="n">warning</span><span class="p">(</span><span class="err">&#39;</span><span class="n">Same</span> <span class="n">label</span> <span class="n">found</span> <span class="k">for</span> <span class="n">multiple</span> <span class="n">electrodes</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC299" class="line">            <span class="n">end</span></span>
<span id="LC300" class="line">        <span class="n">end</span></span>
<span id="LC301" class="line"></span>
<span id="LC302" class="line">        <span class="k">if</span> <span class="p">(</span><span class="n">length</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span></span>
<span id="LC303" class="line">            <span class="n">noteNameNew</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">noteNameNew</span><span class="p">,</span> <span class="n">chLabel1</span><span class="p">,</span><span class="n">csv_all_electrodes</span><span class="p">{</span><span class="n">iChanMatch1</span><span class="p">,</span><span class="mi">1</span><span class="p">});</span> <span class="o">%</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">iChanMatch1</span><span class="p">,</span><span class="mi">2</span><span class="p">}</span></span>
<span id="LC304" class="line">        <span class="n">end</span></span>
<span id="LC305" class="line">        <span class="k">if</span> <span class="p">(</span><span class="n">length</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span></span>
<span id="LC306" class="line">            <span class="n">noteNameNew</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">noteNameNew</span><span class="p">,</span> <span class="n">chLabel2</span><span class="p">,</span><span class="n">csv_all_electrodes</span><span class="p">{</span><span class="n">iChanMatch2</span><span class="p">,</span><span class="mi">1</span><span class="p">});</span> <span class="o">%</span> <span class="n">chMatchLog</span><span class="p">{</span><span class="n">iChanMatch2</span><span class="p">,</span><span class="mi">2</span><span class="p">}</span>            </span>
<span id="LC307" class="line">        <span class="n">end</span></span>
<span id="LC308" class="line">        <span class="k">if</span> <span class="p">(</span><span class="n">length</span><span class="p">(</span><span class="n">iChanMatch1</span><span class="p">)</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span> <span class="o">&amp;&amp;</span>  <span class="p">(</span><span class="n">length</span><span class="p">(</span><span class="n">iChanMatch2</span><span class="p">)</span> <span class="o">==</span> <span class="mi">1</span><span class="p">)</span></span>
<span id="LC309" class="line">        <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">trials</span><span class="p">.</span><span class="n">events</span><span class="p">(</span><span class="n">iEvt</span><span class="p">).</span><span class="n">type</span> <span class="o">=</span> <span class="n">noteNameNew</span><span class="p">;</span></span>
<span id="LC310" class="line">        <span class="n">end</span></span>
<span id="LC311" class="line">    <span class="n">end</span>    </span>
<span id="LC312" class="line">    <span class="n">csv_struct</span><span class="p">.</span><span class="n">csv_labels</span> <span class="o">=</span> <span class="n">csv_all_electrodes</span><span class="p">(</span><span class="o">:</span><span class="p">,</span><span class="mi">1</span><span class="p">);</span></span>
<span id="LC313" class="line">    <span class="n">SpmMat</span><span class="p">.</span><span class="n">D</span><span class="p">.</span><span class="n">other</span> <span class="o">=</span> <span class="n">csv_struct</span><span class="p">;</span> <span class="o">%</span> <span class="n">Add</span> <span class="n">an</span> <span class="n">extra</span> <span class="n">field</span> <span class="n">to</span> <span class="n">the</span> <span class="p">.</span><span class="n">mat</span> <span class="n">so</span> <span class="n">we</span> <span class="n">have</span> <span class="n">a</span> <span class="n">listing</span> <span class="n">of</span> <span class="n">all</span> <span class="n">channel</span> <span class="n">labels</span> <span class="k">in</span> <span class="n">the</span> <span class="n">csv</span><span class="p">.</span></span>
<span id="LC314" class="line">    <span class="o">%</span> <span class="n">Update</span> <span class="n">existing</span> <span class="p">.</span><span class="n">mat</span> <span class="n">file</span></span>
<span id="LC315" class="line">    <span class="n">save</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">,</span> <span class="err">&#39;</span><span class="o">-</span><span class="k">struct</span><span class="err">&#39;</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">SpmMat</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC316" class="line">    <span class="n">save</span><span class="p">(</span><span class="n">spm_eeg_load</span><span class="p">(</span><span class="n">SpmFile</span><span class="p">));</span> <span class="o">%</span> <span class="n">SPM</span><span class="err">&#39;</span><span class="n">s</span> <span class="n">save</span> <span class="n">to</span> <span class="n">create</span> <span class="n">a</span> <span class="n">valid</span> <span class="n">SPM</span> <span class="n">object</span></span>
<span id="LC317" class="line"><span class="n">end</span></span>
<span id="LC318" class="line"></span>
<span id="LC319" class="line"><span class="n">end</span></span>
<span id="LC320" class="line"></span>
<span id="LC321" class="line"></span>
<span id="LC322" class="line"><span class="o">%%</span> <span class="n">Read</span> <span class="n">name</span> <span class="n">from</span> <span class="n">file</span></span>
<span id="LC323" class="line"><span class="n">function</span> <span class="n">Name</span> <span class="o">=</span> <span class="n">readName</span><span class="p">(</span><span class="n">filename</span><span class="p">)</span></span>
<span id="LC324" class="line">    <span class="n">fid</span> <span class="o">=</span> <span class="n">fopen</span><span class="p">(</span><span class="n">filename</span><span class="p">);</span></span>
<span id="LC325" class="line">    <span class="n">i1</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span></span>
<span id="LC326" class="line">    <span class="k">while</span> <span class="mi">1</span></span>
<span id="LC327" class="line">        <span class="n">tmp</span> <span class="o">=</span> <span class="n">fgetl</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC328" class="line">        <span class="k">if</span> <span class="n">isequal</span><span class="p">(</span><span class="n">tmp</span><span class="p">,</span> <span class="o">-</span><span class="mi">1</span><span class="p">)</span> <span class="o">||</span> <span class="n">isempty</span><span class="p">(</span><span class="n">tmp</span><span class="p">)</span></span>
<span id="LC329" class="line">            <span class="k">break</span><span class="p">;</span></span>
<span id="LC330" class="line">        <span class="n">end</span></span>
<span id="LC331" class="line">        <span class="n">Name</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="o">=</span> <span class="n">tmp</span><span class="p">;</span></span>
<span id="LC332" class="line">        <span class="n">i1</span> <span class="o">=</span> <span class="n">i1</span> <span class="o">+</span> <span class="mi">1</span><span class="p">;</span></span>
<span id="LC333" class="line">    <span class="n">end</span></span>
<span id="LC334" class="line">    <span class="n">fclose</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC335" class="line"><span class="n">end</span></span>
<span id="LC336" class="line"></span>
<span id="LC337" class="line"></span>
<span id="LC338" class="line"><span class="o">%%</span> <span class="n">Read</span> <span class="n">name</span> <span class="n">and</span> <span class="n">position</span> <span class="n">from</span> <span class="p">.</span><span class="n">csv</span></span>
<span id="LC339" class="line"><span class="n">function</span> <span class="p">[</span><span class="n">Name</span><span class="p">,</span> <span class="nf">Pos</span><span class="p">]</span> <span class="o">=</span> <span class="n">readCsv</span><span class="p">(</span><span class="n">filename</span><span class="p">)</span></span>
<span id="LC340" class="line">    <span class="n">Name</span> <span class="o">=</span> <span class="p">{};</span></span>
<span id="LC341" class="line">    <span class="n">Pos</span> <span class="o">=</span> <span class="p">[];</span></span>
<span id="LC342" class="line">    <span class="o">%</span> <span class="n">Open</span> <span class="n">file</span></span>
<span id="LC343" class="line">    <span class="n">fid</span> <span class="o">=</span> <span class="n">fopen</span><span class="p">(</span><span class="n">filename</span><span class="p">);</span></span>
<span id="LC344" class="line">    <span class="n">i1</span> <span class="o">=</span> <span class="mi">1</span><span class="p">;</span></span>
<span id="LC345" class="line">    <span class="o">%</span> <span class="n">Skip</span> <span class="n">two</span> <span class="n">header</span> <span class="n">lines</span></span>
<span id="LC346" class="line">    <span class="n">fgetl</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC347" class="line">    <span class="n">fgetl</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC348" class="line">    <span class="o">%</span> <span class="n">Read</span> <span class="n">column</span> <span class="n">names</span> <span class="p">(</span><span class="n">tab</span><span class="o">-</span><span class="n">separated</span> <span class="n">file</span><span class="p">)</span></span>
<span id="LC349" class="line">    <span class="n">colNames</span> <span class="o">=</span> <span class="n">strsplit</span><span class="p">(</span><span class="n">fgetl</span><span class="p">(</span><span class="n">fid</span><span class="p">),</span> <span class="sc">&#39;\t&#39;</span><span class="p">);</span></span>
<span id="LC350" class="line">    <span class="n">iColName</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">colNames</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">contact</span><span class="err">&#39;</span><span class="p">));</span></span>
<span id="LC351" class="line">    <span class="n">iColPos</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmpi</span><span class="p">(</span><span class="n">colNames</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">T1pre</span> <span class="n">Scanner</span> <span class="n">Based</span><span class="err">&#39;</span><span class="p">));</span></span>
<span id="LC352" class="line">    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iColName</span><span class="p">)</span> <span class="o">||</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iColPos</span><span class="p">)</span></span>
<span id="LC353" class="line">        <span class="n">disp</span><span class="p">(</span><span class="err">&#39;</span><span class="n">ImaGIN</span><span class="o">&gt;</span> <span class="n">ERROR</span><span class="o">:</span> <span class="n">Invalid</span> <span class="n">CSV</span> <span class="n">file</span><span class="o">:</span> <span class="n">No</span> <span class="n">columns</span> <span class="s">&quot;T1pre Scanner Based&quot;</span> <span class="n">or</span> <span class="s">&quot;contact&quot;</span><span class="p">.</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC354" class="line">        <span class="k">return</span><span class="p">;</span></span>
<span id="LC355" class="line">    <span class="n">end</span></span>
<span id="LC356" class="line">    <span class="o">%</span> <span class="n">Read</span> <span class="n">contact</span> <span class="n">by</span> <span class="n">contact</span></span>
<span id="LC357" class="line">    <span class="k">while</span> <span class="mi">1</span></span>
<span id="LC358" class="line">        <span class="o">%</span> <span class="n">Read</span> <span class="n">line</span></span>
<span id="LC359" class="line">        <span class="n">tmp</span> <span class="o">=</span> <span class="n">fgetl</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC360" class="line">        <span class="k">if</span> <span class="n">isequal</span><span class="p">(</span><span class="n">tmp</span><span class="p">,</span> <span class="o">-</span><span class="mi">1</span><span class="p">)</span> <span class="o">||</span> <span class="n">isempty</span><span class="p">(</span><span class="n">tmp</span><span class="p">)</span></span>
<span id="LC361" class="line">            <span class="k">break</span><span class="p">;</span></span>
<span id="LC362" class="line">        <span class="n">end</span></span>
<span id="LC363" class="line">        <span class="o">%</span> <span class="n">Split</span> <span class="n">by</span> <span class="n">columns</span> <span class="p">(</span><span class="n">tab</span><span class="o">-</span><span class="n">separated</span> <span class="n">file</span><span class="p">)</span></span>
<span id="LC364" class="line">        <span class="n">tmp</span> <span class="o">=</span> <span class="n">strsplit</span><span class="p">(</span><span class="n">tmp</span><span class="p">,</span> <span class="sc">&#39;\t&#39;</span><span class="p">);</span></span>
<span id="LC365" class="line">        <span class="n">Name</span><span class="p">{</span><span class="n">i1</span><span class="p">}</span> <span class="o">=</span> <span class="n">tmp</span><span class="p">{</span><span class="n">iColName</span><span class="p">};</span></span>
<span id="LC366" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">str2num</span><span class="p">(</span><span class="n">tmp</span><span class="p">{</span><span class="n">iColPos</span><span class="p">}))</span></span>
<span id="LC367" class="line">            <span class="n">Pos</span><span class="p">(</span><span class="n">i1</span><span class="p">,</span><span class="mi">1</span><span class="o">:</span><span class="mi">3</span><span class="p">)</span> <span class="o">=</span> <span class="n">str2num</span><span class="p">(</span><span class="n">tmp</span><span class="p">{</span><span class="n">iColPos</span><span class="p">});</span></span>
<span id="LC368" class="line">        <span class="k">else</span></span>
<span id="LC369" class="line">            <span class="n">Pos</span><span class="p">(</span><span class="n">i1</span><span class="p">,</span><span class="mi">1</span><span class="o">:</span><span class="mi">3</span><span class="p">)</span> <span class="o">=</span> <span class="n">NaN</span><span class="p">;</span></span>
<span id="LC370" class="line">        <span class="n">end</span></span>
<span id="LC371" class="line">        <span class="n">i1</span> <span class="o">=</span> <span class="n">i1</span> <span class="o">+</span> <span class="mi">1</span><span class="p">;</span></span>
<span id="LC372" class="line">    <span class="n">end</span></span>
<span id="LC373" class="line">    <span class="n">fclose</span><span class="p">(</span><span class="n">fid</span><span class="p">);</span></span>
<span id="LC374" class="line"><span class="n">end</span></span>
<span id="LC375" class="line"></span>
<span id="LC376" class="line"></span>
<span id="LC377" class="line"><span class="o">%%</span> <span class="n">Find</span> <span class="n">channel</span> <span class="n">name</span> <span class="k">in</span> <span class="n">a</span> <span class="n">list</span></span>
<span id="LC378" class="line"><span class="n">function</span> <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">findChannel</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="n">List</span><span class="p">,</span> <span class="n">caseType</span><span class="p">)</span></span>
<span id="LC379" class="line">    <span class="o">%</span> <span class="n">Test</span> <span class="n">three</span> <span class="n">different</span> <span class="k">case</span> <span class="n">versions</span></span>
<span id="LC380" class="line">    <span class="o">%</span> <span class="n">The</span> <span class="n">goal</span> <span class="n">is</span> <span class="n">to</span> <span class="n">be</span> <span class="n">able</span> <span class="n">to</span> <span class="n">handle</span> <span class="n">a</span> <span class="n">List</span> <span class="n">containing</span> <span class="n">both</span> <span class="n">electrodes</span> <span class="n">L</span><span class="err">&#39;</span> <span class="p">(</span><span class="s">&quot;Lp&quot;</span><span class="p">)</span> <span class="n">and</span> <span class="n">LP</span> <span class="p">(</span><span class="s">&quot;LP&quot;</span><span class="p">)</span></span>
<span id="LC381" class="line">    <span class="o">%</span> <span class="n">The</span> <span class="n">List</span> <span class="n">provided</span> <span class="n">here</span> <span class="n">comes</span> <span class="n">from</span> <span class="n">a</span> <span class="p">.</span><span class="n">csv</span> <span class="n">with</span> <span class="n">this</span> <span class="n">convention</span><span class="p">:</span> <span class="n">only</span> <span class="n">capitals</span> <span class="n">except</span> <span class="k">for</span> <span class="n">the</span> <span class="s">&quot;p&quot;</span> <span class="k">for</span> <span class="s">&quot;&#39;&quot;</span></span>
<span id="LC382" class="line">    <span class="k">if</span> <span class="p">(</span><span class="n">nargin</span> <span class="o">&lt;</span> <span class="mi">3</span><span class="p">)</span> <span class="o">||</span> <span class="n">isempty</span><span class="p">(</span><span class="n">caseType</span><span class="p">)</span></span>
<span id="LC383" class="line">        <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">findChannel</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="n">List</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">no_change</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC384" class="line">        <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC385" class="line">            <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">findChannel</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="n">List</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">all_upper</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC386" class="line">        <span class="n">end</span></span>
<span id="LC387" class="line">        <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC388" class="line">            <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">findChannel</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="n">List</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">upper_except_p</span><span class="err">&#39;</span><span class="p">);</span></span>
<span id="LC389" class="line">        <span class="n">end</span></span>
<span id="LC390" class="line">        <span class="k">return</span><span class="p">;</span></span>
<span id="LC391" class="line">    <span class="n">end</span></span>
<span id="LC392" class="line">    </span>
<span id="LC393" class="line">    <span class="o">%</span> <span class="n">Remove</span> <span class="n">spaces</span></span>
<span id="LC394" class="line">    <span class="n">Label</span><span class="p">(</span><span class="n">Label</span> <span class="o">==</span> <span class="sc">&#39; &#39;</span><span class="p">)</span> <span class="o">=</span> <span class="p">[];</span></span>
<span id="LC395" class="line">    <span class="o">%</span> <span class="n">Remove</span> <span class="n">spaces</span></span>
<span id="LC396" class="line">    <span class="n">List</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">List</span><span class="p">,</span><span class="sc">&#39;-&#39;</span><span class="p">,</span><span class="err">&#39;&#39;</span><span class="p">);</span>    </span>
<span id="LC397" class="line">    <span class="o">%</span> <span class="n">Switch</span> <span class="k">case</span> <span class="n">type</span></span>
<span id="LC398" class="line">    <span class="k">switch</span> <span class="p">(</span><span class="n">caseType</span><span class="p">)</span></span>
<span id="LC399" class="line">        <span class="k">case</span> <span class="err">&#39;</span><span class="n">no_change</span><span class="err">&#39;</span></span>
<span id="LC400" class="line">            <span class="o">%</span> <span class="n">Nothing</span> <span class="n">to</span> <span class="n">change</span></span>
<span id="LC401" class="line">        <span class="k">case</span> <span class="err">&#39;</span><span class="n">all_upper</span><span class="err">&#39;</span></span>
<span id="LC402" class="line">            <span class="n">Label</span> <span class="o">=</span> <span class="n">upper</span><span class="p">(</span><span class="n">Label</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC403" class="line">        <span class="k">case</span> <span class="err">&#39;</span><span class="n">upper_except_p</span><span class="err">&#39;</span></span>
<span id="LC404" class="line">            <span class="n">iP</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">Label</span> <span class="o">==</span> <span class="sc">&#39;p&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC405" class="line">            <span class="n">Label</span> <span class="o">=</span> <span class="n">upper</span><span class="p">(</span><span class="n">Label</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC406" class="line">            <span class="k">if</span> <span class="n">any</span><span class="p">(</span><span class="n">iP</span> <span class="o">&gt;=</span> <span class="mi">2</span><span class="p">)</span></span>
<span id="LC407" class="line">                <span class="n">Label</span><span class="p">(</span><span class="n">iP</span><span class="p">(</span><span class="n">end</span><span class="p">))</span> <span class="o">=</span> <span class="sc">&#39;p&#39;</span><span class="err">;</span></span>
<span id="LC408" class="line">            <span class="n">end</span></span>
<span id="LC409" class="line">    <span class="n">end</span></span>
<span id="LC410" class="line">    <span class="o">%</span> <span class="n">Replacing</span> <span class="err">&#39;</span> <span class="n">with</span> <span class="n">p</span></span>
<span id="LC411" class="line">    <span class="n">Label</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="err">&#39;&#39;&#39;&#39;</span><span class="p">,</span> <span class="sc">&#39;p&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC412" class="line">    <span class="n">List</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">List</span><span class="p">,</span> <span class="err">&#39;&#39;&#39;&#39;</span><span class="p">,</span> <span class="sc">&#39;p&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC413" class="line">    <span class="o">%</span> <span class="n">Replacing</span> <span class="p">,</span> <span class="n">with</span> <span class="n">p</span></span>
<span id="LC414" class="line">    <span class="n">Label</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="sc">&#39;,&#39;</span><span class="p">,</span> <span class="sc">&#39;p&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC415" class="line">    <span class="n">List</span> <span class="o">=</span> <span class="n">strrep</span><span class="p">(</span><span class="n">List</span><span class="p">,</span> <span class="sc">&#39;,&#39;</span><span class="p">,</span> <span class="sc">&#39;p&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC416" class="line">    <span class="o">%</span> <span class="n">Look</span> <span class="k">for</span> <span class="n">channel</span> <span class="k">in</span> <span class="n">position</span> <span class="n">file</span></span>
<span id="LC417" class="line">    <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmp</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="n">List</span><span class="p">))</span><span class="err">;</span></span>
<span id="LC418" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC419" class="line">        <span class="k">return</span><span class="err">;</span></span>
<span id="LC420" class="line">    <span class="n">end</span>    </span>
<span id="LC421" class="line">    <span class="o">%</span> <span class="n">Channel</span> <span class="n">not</span> <span class="n">found</span><span class="p">:</span> <span class="n">try</span> <span class="n">adding</span> <span class="n">missing</span> <span class="mi">0</span> <span class="p">(</span><span class="n">A1</span><span class="o">=&gt;</span><span class="n">A01</span><span class="p">,</span> <span class="n">T11</span><span class="o">=&gt;</span><span class="n">T101</span><span class="p">)</span></span>
<span id="LC422" class="line">    <span class="o">%</span> <span class="n">Find</span> <span class="n">the</span> <span class="n">last</span> <span class="n">letter</span> <span class="k">in</span> <span class="n">the</span> <span class="n">name</span></span>
<span id="LC423" class="line">    <span class="n">iLastLetter</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="o">~</span><span class="n">ismember</span><span class="p">(</span><span class="n">Label</span><span class="p">,</span> <span class="err">&#39;</span><span class="mo">01234567</span><span class="mi">89</span><span class="err">&#39;</span><span class="p">),</span> <span class="mi">1</span><span class="p">,</span> <span class="err">&#39;</span><span class="n">last</span><span class="err">&#39;</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC424" class="line">    <span class="o">%</span> <span class="n">If</span> <span class="n">there</span> <span class="n">is</span> <span class="n">not</span> <span class="n">electrode</span> <span class="n">label</span> <span class="n">or</span> <span class="n">no</span> <span class="n">index</span><span class="p">:</span> <span class="n">wrong</span> <span class="n">naming</span></span>
<span id="LC425" class="line">    <span class="k">if</span> <span class="n">isempty</span><span class="p">(</span><span class="n">iLastLetter</span><span class="p">)</span> <span class="o">||</span> <span class="p">(</span><span class="n">iLastLetter</span> <span class="o">==</span> <span class="n">length</span><span class="p">(</span><span class="n">Label</span><span class="p">))</span></span>
<span id="LC426" class="line">        <span class="k">return</span><span class="err">;</span></span>
<span id="LC427" class="line">    <span class="n">end</span></span>
<span id="LC428" class="line">    <span class="o">%</span> <span class="n">Find</span> <span class="n">channel</span> <span class="n">index</span></span>
<span id="LC429" class="line">    <span class="n">chLabel</span> <span class="o">=</span> <span class="n">Label</span><span class="p">(</span><span class="mi">1</span><span class="p">:</span><span class="n">iLastLetter</span><span class="p">)</span><span class="err">;</span></span>
<span id="LC430" class="line">    <span class="n">chInd</span> <span class="o">=</span> <span class="n">str2num</span><span class="p">(</span><span class="n">Label</span><span class="p">(</span><span class="n">iLastLetter</span><span class="o">+</span><span class="mi">1</span><span class="p">:</span><span class="n">end</span><span class="p">));</span></span>
<span id="LC431" class="line">    <span class="o">%</span> <span class="n">If</span> <span class="n">label</span> <span class="n">is</span> <span class="o">&gt;</span> <span class="mi">100</span><span class="o">:</span> <span class="n">Add</span> <span class="mi">1</span> <span class="n">or</span> <span class="mi">2</span> <span class="n">at</span> <span class="n">the</span> <span class="n">end</span> <span class="n">of</span> <span class="n">the</span> <span class="n">label</span></span>
<span id="LC432" class="line">    <span class="k">if</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&gt;</span> <span class="mi">220</span><span class="p">)</span></span>
<span id="LC433" class="line">        <span class="k">return</span><span class="p">;</span></span>
<span id="LC434" class="line">    <span class="n">elseif</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&gt;</span> <span class="mi">200</span><span class="p">)</span></span>
<span id="LC435" class="line">        <span class="n">chLabel</span> <span class="o">=</span> <span class="p">[</span><span class="n">chLabel</span> <span class="sc">&#39;2&#39;</span><span class="p">];</span></span>
<span id="LC436" class="line">        <span class="n">chInd</span> <span class="o">=</span> <span class="n">chInd</span> <span class="o">-</span> <span class="mi">200</span><span class="p">;</span></span>
<span id="LC437" class="line">    <span class="n">elseif</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&gt;</span> <span class="mi">100</span><span class="p">)</span></span>
<span id="LC438" class="line">        <span class="n">chLabel</span> <span class="o">=</span> <span class="p">[</span><span class="n">chLabel</span> <span class="sc">&#39;1&#39;</span><span class="p">];</span></span>
<span id="LC439" class="line">        <span class="n">chInd</span> <span class="o">=</span> <span class="n">chInd</span> <span class="o">-</span> <span class="mi">100</span><span class="p">;</span></span>
<span id="LC440" class="line">    <span class="n">elseif</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&gt;</span> <span class="mi">20</span><span class="p">)</span></span>
<span id="LC441" class="line">        <span class="n">chLabel</span> <span class="o">=</span> <span class="p">[</span><span class="n">chLabel</span> <span class="sc">&#39;2&#39;</span><span class="p">];</span></span>
<span id="LC442" class="line">        <span class="n">chInd</span> <span class="o">=</span> <span class="n">chInd</span> <span class="o">-</span> <span class="mi">20</span><span class="p">;</span></span>
<span id="LC443" class="line">    <span class="n">end</span></span>
<span id="LC444" class="line">    </span>
<span id="LC445" class="line">    <span class="o">%</span> <span class="n">Try</span> <span class="n">with</span> <span class="n">a</span> <span class="n">zero</span></span>
<span id="LC446" class="line">    <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmp</span><span class="p">(</span><span class="n">sprintf</span><span class="p">(</span><span class="err">&#39;</span><span class="o">%</span><span class="n">s</span><span class="o">%</span><span class="mo">02</span><span class="n">d</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chLabel</span><span class="p">,</span> <span class="n">chInd</span><span class="p">),</span> <span class="n">List</span><span class="p">));</span></span>
<span id="LC447" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC448" class="line">        <span class="k">return</span><span class="p">;</span></span>
<span id="LC449" class="line">    <span class="n">end</span></span>
<span id="LC450" class="line">    <span class="o">%</span> <span class="n">Try</span> <span class="n">without</span> <span class="n">the</span> <span class="n">zero</span></span>
<span id="LC451" class="line">    <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmp</span><span class="p">(</span><span class="n">sprintf</span><span class="p">(</span><span class="err">&#39;</span><span class="o">%</span><span class="n">s</span><span class="o">%</span><span class="n">d</span><span class="err">&#39;</span><span class="p">,</span> <span class="n">chLabel</span><span class="p">,</span> <span class="n">chInd</span><span class="p">),</span> <span class="n">List</span><span class="p">));</span></span>
<span id="LC452" class="line">    <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC453" class="line">        <span class="k">return</span><span class="p">;</span></span>
<span id="LC454" class="line">    <span class="n">end</span></span>
<span id="LC455" class="line"></span>
<span id="LC456" class="line">    <span class="o">%</span> <span class="n">For</span> <span class="n">channel</span> <span class="n">indices</span> <span class="n">between</span> <span class="mi">11</span> <span class="n">and</span> <span class="mi">19</span><span class="p">,</span> <span class="n">try</span> <span class="n">to</span> <span class="n">interpret</span> <span class="n">them</span> <span class="n">as</span> <span class="k">if</span> <span class="n">electrode</span> <span class="n">name</span> <span class="n">was</span> <span class="n">ending</span> <span class="n">with</span> <span class="n">a</span> <span class="mi">1</span></span>
<span id="LC457" class="line">    <span class="k">if</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&gt;=</span> <span class="mi">11</span><span class="p">)</span> <span class="o">&amp;&amp;</span> <span class="p">(</span><span class="n">chInd</span> <span class="o">&lt;=</span> <span class="mi">19</span><span class="p">)</span></span>
<span id="LC458" class="line">        <span class="n">iChanPos</span> <span class="o">=</span> <span class="n">find</span><span class="p">(</span><span class="n">strcmp</span><span class="p">(</span><span class="n">sprintf</span><span class="p">(</span><span class="err">&#39;</span><span class="o">%</span><span class="n">s</span><span class="o">%</span><span class="mo">02</span><span class="n">d</span><span class="err">&#39;</span><span class="p">,</span> <span class="p">[</span><span class="n">chLabel</span> <span class="sc">&#39;1&#39;</span><span class="p">],</span> <span class="n">chInd</span> <span class="o">-</span> <span class="mi">10</span><span class="p">),</span> <span class="n">List</span><span class="p">));</span></span>
<span id="LC459" class="line">        <span class="k">if</span> <span class="o">~</span><span class="n">isempty</span><span class="p">(</span><span class="n">iChanPos</span><span class="p">)</span></span>
<span id="LC460" class="line">            <span class="k">return</span><span class="p">;</span></span>
<span id="LC461" class="line">        <span class="n">end</span></span>
<span id="LC462" class="line">    <span class="n">end</span></span>
<span id="LC463" class="line"><span class="n">end</span></span>
<span id="LC464" class="line"></span></code></pre>

</div>
</div>


</article>
</div>

</div>

</div>
</div>
</div>
</div>
</div>



</body>
</html>

