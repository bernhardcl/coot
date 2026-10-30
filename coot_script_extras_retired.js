// Based on Sharex github pages: https://getsharex.com/downloads/

$(document).ready(function () {
      GetLatestReleaseInfo();
   });

function GetLatestReleaseInfo() {

   $.getJSON("https://api.github.com/repos/bernhardcl/coot/releases/latest").done(function (release) {
         var asset = release.assets[0];
         var downloadCount = 0;
         var release_name;
         for (var i = 0; i < release.assets.length; i++) {
            downloadCount += release.assets[i].download_count;
         }
         var oneHour = 60 * 60 * 1000;
         var oneDay = 24 * oneHour;
         var dateDiff = new Date() - new Date(asset.updated_at);
         var timeAgo;
         if (dateDiff < oneDay)
            {
               timeAgo = (dateDiff / oneHour).toFixed(1) + " hours ago";
            }
         else
            {
               timeAgo = (dateDiff / oneDay).toFixed(1) + " days ago";
            }

         var releaseInfo = release.name + " was updated " + timeAgo + " and downloaded " + downloadCount.toLocaleString() + " times.";
         $(".coot-download-href").attr("href", asset.browser_download_url);
         $(".coot-download-action").attr("action", asset.browser_download_url);
         $(".coot-download-info").attr("value", "Download " + release.name);
         $(".release-info").text(releaseInfo);
         $(".version-info").text(release.name);
      });
}


function GetReleases(repo, preReleaseFilter, versionFilter, showDownloadCounts) {
   // preReleaseFilter: "releases" (only stable), "prereleases" (only pre-releases), "all" (both)
   // versionFilter: "all", "0.9", "1.x" (or "latest" for non-0.9)
   // showDownloadCounts: true/false (default true)

   if (typeof showDownloadCounts === 'undefined') {
      showDownloadCounts = true;
   }

   $.getJSON("https://api.github.com/repos/" + repo + "/releases").done(function(json) {
      var totalDownloadCount = 0;

      for (var i = 0; i < json.length; i++) {
         var release = json[i];

         if (release.assets.length === 0) {
            continue;
         }

         // Pre-release filtering
         var includeByPreRelease = false;
         if (preReleaseFilter === "all") {
            includeByPreRelease = true;
         } else if (preReleaseFilter === "releases") {
            includeByPreRelease = !release.prerelease;
         } else if (preReleaseFilter === "prereleases") {
            includeByPreRelease = release.prerelease;
         }

         if (!includeByPreRelease) {
            continue;
         }

         // Version filtering
         var includeByVersion = false;
         // Extract major version (first number after common prefixes)
         var versionMatch = release.tag_name.match(/(\d+)\.\d+/) ||
                            release.name.match(/(\d+)\.\d+/);
         var majorVersion = versionMatch ? versionMatch[1] : null;
         if (versionFilter === "all") {
            includeByVersion = true;
         } else if (versionFilter === "0.9") {
            includeByVersion = (majorVersion === "0");
         } else if (versionFilter === "1.x" || versionFilter === "latest") {
            includeByVersion = (majorVersion === "1");
         }

         if (!includeByVersion) {
            continue;
         }

         var asset = release.assets[0];
         var fileSize = Math.round(asset.size / (1024*1024));
         var downloadCount = 0;

         for (var i2 = 0; i2 < release.assets.length; i2++) {
            downloadCount += release.assets[i2].download_count;
         }

         var md5sum_url = "dummy";
         for (var j = 0; j < release.assets.length; j++) {
            var asset2 = release.assets[j];
            var md5sum = asset.name + ".md5sum";
            if (asset2.name == md5sum) {
               md5sum_url = asset2.browser_download_url;
            }
         }

         totalDownloadCount += downloadCount;

         var d = new Date(asset.updated_at);

         var row = $("<tr>")
            .append($("<td>")
               .append($("<a>")
                  .attr("href", release.html_url)
                  .text(release.name)
               )
            )
            .append($("<td>")
               .append($("<a>")
                  .attr("href", asset.browser_download_url)
                  .text("download installer")
               )
            )
            .append($("<td>")
               .append($("<a>")
                  .attr("href", md5sum_url)
                  .text("get md5sum")
               )
            )
            .append($("<td>")
               .append($("<a>")
                  .text(fileSize.toLocaleString() + " MB")
               )
            )
            .append($("<td>")
               .text(d.toLocaleString())
            );

         // Conditionally add download count column
         if (showDownloadCounts) {
            row.append($("<td>")
               .text(downloadCount.toLocaleString())
            );
         }

         $(".table-downloads tbody").append(row);
      }

      $(".table-downloads").show();
   });
}
