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
                           
         $(".release-info").text(releaseInfo);
         $(".version-info").text(release.name);               

      });
}
      
function GetReleases(repo, pre_release) {
   $.getJSON("https://api.github.com/repos/" + repo + "/releases").done(function(json) {
         var totalDownloadCount = 0;
         for (var i = 0; i < json.length; i++) {
            var release = json[i];
            if (release.assets.length === 0) {
               continue;
            }
            var asset = release.assets[0];
            var fileSize = Math.round(asset.size / (1024*1024));
            var downloadCount = 0;
            for (var i2 = 0; i2 < release.assets.length; i2++) {
               downloadCount += release.assets[i2].download_count;
            }
            totalDownloadCount += downloadCount;
            // to avoid an extra library (moment.js)
            var d = new Date(asset.updated_at);
            if (pre_release == release.prerelease) {
               $(".table-downloads tbody")
                  .append($("<tr>")
                          .append($("<td>")
                                  .append($("<a>")
                                          .attr("href", asset.browser_download_url)
                                          .text(release.name)
                                          )
                                  )
                          .append($("<td>")
                                  .append($("<a>")
                                          //                                .attr("href", asset.browser_download_url)
                                          .text(fileSize.toLocaleString() + " MB")
                                          )
                                  )
                          //                        .append($("<td>")
                          //                            .text(downloadCount.toLocaleString())
                          //                            )
                          .append($("<td>")
                                  .text(d.toLocaleString())
                                  //                            .text(moment(asset.updated_at).format("YYYY-MM-DD HH:mm"))
                                  )
                          );
            }
         }

         //         if (totalDownloadCount > 0)
         //  {
         //      $(".total-downloads").text(" (" + totalDownloadCount.toLocaleString() + " downloads)");
         //   }

         $(".table-downloads").show();
         //         console.log($(".table-downloads").prop('outerHTML'));
      });
}

function change_myselect(sel) {
   // clean the table
   $(".table-downloads tbody").empty();

   $(document).ready(function() {
         repo = "bernhardcl/coot";
         if (sel == "all") {
            GetReleases(repo, false);
            GetReleases(repo, true);
         } else {
            var isTrueSet = (sel == 'true');
            GetReleases(repo, isTrueSet);
         }
      });
}
