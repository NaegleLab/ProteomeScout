<?php
header("HTTP/1.1 503 Service Temporarily Unavailable");
header("Status: 503 Service Temporarily Unavailable");
header("Retry-After: 3600");
?>
<html>
    <head>
        <title>PTMScout Temporarily Unavailable</title>
    </head>
    <body>
        <h1>PTMScout 2.0 is down for maintenance.</h1>
        <p>We are hard at work improving PTMScout. Please try your request again at a later time.</p>
    </body>
</html>
