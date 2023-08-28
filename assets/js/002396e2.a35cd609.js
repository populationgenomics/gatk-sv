"use strict";(self.webpackChunkGATK_SV=self.webpackChunkGATK_SV||[]).push([[4123],{3905:(e,t,n)=>{n.d(t,{Zo:()=>s,kt:()=>g});var i=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function a(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);t&&(i=i.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,i)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?a(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,i,r=function(e,t){if(null==e)return{};var n,i,r={},a=Object.keys(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var c=i.createContext({}),d=function(e){var t=i.useContext(c),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},s=function(e){var t=d(e.components);return i.createElement(c.Provider,{value:t},e.children)},m="mdxType",p={inlineCode:"code",wrapper:function(e){var t=e.children;return i.createElement(i.Fragment,{},t)}},u=i.forwardRef((function(e,t){var n=e.components,r=e.mdxType,a=e.originalType,c=e.parentName,s=l(e,["components","mdxType","originalType","parentName"]),m=d(n),u=r,g=m["".concat(c,".").concat(u)]||m[u]||p[u]||a;return n?i.createElement(g,o(o({ref:t},s),{},{components:n})):i.createElement(g,o({ref:t},s))}));function g(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var a=n.length,o=new Array(a);o[0]=u;var l={};for(var c in t)hasOwnProperty.call(t,c)&&(l[c]=t[c]);l.originalType=e,l[m]="string"==typeof e?e:r,o[1]=l;for(var d=2;d<a;d++)o[d]=n[d];return i.createElement.apply(null,o)}return i.createElement.apply(null,n)}u.displayName="MDXCreateElement"},5718:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>o,default:()=>p,frontMatter:()=>a,metadata:()=>l,toc:()=>d});var i=n(7462),r=(n(7294),n(3905));const a={title:"Incremental Publishing",description:"Incremental Publishing Strategy",sidebar_position:4},o=void 0,l={unversionedId:"advanced/docker/deploy/incremental",id:"advanced/docker/deploy/incremental",title:"Incremental Publishing",description:"Incremental Publishing Strategy",source:"@site/docs/advanced/docker/deploy/incremental.md",sourceDirName:"advanced/docker/deploy",slug:"/advanced/docker/deploy/incremental",permalink:"/gatk-sv/docs/advanced/docker/deploy/incremental",draft:!1,editUrl:"https://github.com/broadinstitute/gatk-sv/tree/master/website/docs/advanced/docker/deploy/incremental.md",tags:[],version:"current",sidebarPosition:4,frontMatter:{title:"Incremental Publishing",description:"Incremental Publishing Strategy",sidebar_position:4},sidebar:"tutorialSidebar",previous:{title:"Manual Deployment",permalink:"/gatk-sv/docs/advanced/docker/deploy/manual"}},c={},d=[{value:"Determining Modified Files",id:"determining-modified-files",level:2},{value:"Identifying Images Requiring Rebuilding from Changed Files",id:"identifying-images-requiring-rebuilding-from-changed-files",level:2}],s={toc:d},m="wrapper";function p(e){let{components:t,...n}=e;return(0,r.kt)(m,(0,i.Z)({},s,n,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("p",null,"The hierarchical and modular organization of GATK-SV Docker\nimages offers a significant advantage: when updating the codebase,\nnot every Docker image is affected, minimizing the impact of changes.\nThis means that not all Docker images need to be rebuilt and\npublished with each pipeline modification. The\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/gatk-sv/blob/main/scripts/docker/build_docker.py"},(0,r.kt)("inlineCode",{parentName:"a"},"build_docker")),"\nscript efficiently tracks these changes and determines which\nDocker images are impacted. Consequently, only the affected Docker\nimages are built, saving both storage space and build time."),(0,r.kt)("p",null,"This incremental and selective building and publishing\nstrategy is particularly beneficial considering the size and\nbuild time of Docker images. By building and publishing\nonly the necessary images, we can save on storage space and\nreduce the overall build time.\nThis page provides a detailed explanation of\nthis incremental and selective approach."),(0,r.kt)("h2",{id:"determining-modified-files"},"Determining Modified Files"),(0,r.kt)("p",null,"The incremental build strategy relies on the determination\nof modified files to identify which Docker images require rebuilding.\nUsing ",(0,r.kt)("inlineCode",{parentName:"p"},"git")," history, the ",(0,r.kt)("inlineCode",{parentName:"p"},"build_docker")," script automatically\ninfers the list of changed files."),(0,r.kt)("p",null,"To achieve this, the script compares two\n",(0,r.kt)("a",{parentName:"p",href:"https://docs.github.com/en/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/about-commits"},(0,r.kt)("inlineCode",{parentName:"a"},"git")," commit SHAs"),": "),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"BASE_SHA"),": the reference commit representing the base branch\n(e.g., ",(0,r.kt)("inlineCode",{parentName:"li"},"broadinstitute/gatk-sv:main"),"), and;"),(0,r.kt)("li",{parentName:"ul"},(0,r.kt)("inlineCode",{parentName:"li"},"HEAD_SHA"),": the target commit representing the latest commit on the feature branch.\n")),(0,r.kt)("p",null,"By analyzing the changes between these commits\nthe script identifies the impacted files and proceeds to\nbuild the corresponding Docker images."),(0,r.kt)("p",null,"During manual runs, the user provides the commit SHAs,\nwhile in automated builds as part of CI/CD,\nthe commit SHAs are determined automatically. "),(0,r.kt)("p",null,"In CI/CD, the commit SHAs are determined as the following example."),(0,r.kt)("mermaid",{value:'%%{init: { \n            \'logLevel\': \'debug\',\n            \'gitGraph\': {\'rotateCommitLabel\': false}, \n            \'themeVariables\': { \'commitLabelFontSize\': \'22px\' } \n         } \n   }%%\ngitGraph\n   commit id: "A"\n   commit id: "B"\n   branch feature\n   checkout feature\n   commit id: "X"\n   checkout main\n   commit id: "C"\n   checkout feature\n   commit id: "Y"\n   checkout main\n   commit id: "D"\n   checkout feature\n   commit id: "Z"\n   checkout main\n   merge feature id: "E"\n   commit id: "F"'}),(0,r.kt)("p",null,"In this example, ",(0,r.kt)("inlineCode",{parentName:"p"},"BASE_SHA=B"),", ",(0,r.kt)("inlineCode",{parentName:"p"},"HEAD_SHA=Z"),", and ",(0,r.kt)("inlineCode",{parentName:"p"},"E")," is the merge commit."),(0,r.kt)("h2",{id:"identifying-images-requiring-rebuilding-from-changed-files"},"Identifying Images Requiring Rebuilding from Changed Files"),(0,r.kt)("p",null,"The ",(0,r.kt)("inlineCode",{parentName:"p"},"build_docker")," script identifies the list of docker images\nthat need to be rebuilt based on two factors. "),(0,r.kt)("ol",null,(0,r.kt)("li",{parentName:"ol"},(0,r.kt)("p",{parentName:"li"},"Directly impacted images are determined by checking the\nlist of files each image depends on. If any of these files have\nchanged, the corresponding image needs rebuilding. ")),(0,r.kt)("li",{parentName:"ol"},(0,r.kt)("p",{parentName:"li"},"Indirectly impacted images are identified based on\nthe hierarchical dependency between images.\nIf a base image is rebuilt, any dependent images built upon\nit also require rebuilding.\n"))),(0,r.kt)("p",null,"This two-step process ensures that all the affected images are correctly\nidentified for rebuilding."),(0,r.kt)("p",null,"A comprehensive mapping of files to their corresponding\nDocker images, specifying which images need to be\nrebuilt when their associated files are updated is given in\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/broadinstitute/gatk-sv/blob/e86d59962146ae1770c535a97c2774d825026957/scripts/docker/build_docker.py#L170-L245"},"this section"),"."))}p.isMDXComponent=!0}}]);