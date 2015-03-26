function highlightEss(i) {
  var text = $(this).text();
  var f = parseFloat(text);
  if(f < 200.0) {
    $(this).append("<i class='icon-asterisk'></i>");
    if(f < 100.0) $(this).addClass('ess-low');
    else $(this).addClass('ess-medium');
  }
}

$(function() {
  $("td.statistic-ess").each(highlightEss);
})
