(function () {
  "use strict";

  document.addEventListener("DOMContentLoaded", function () {
    var configEl = document.getElementById("documenter-copybutton-config");
    if (!configEl) return;

    var config;
    try {
      config = JSON.parse(configEl.textContent);
    } catch (e) {
      console.error("DocumenterCopyButton: failed to parse config", e);
      return;
    }

    var mdUrl = config.mdUrl || "";
    var providers = config.providers || [];
    var prompt = config.prompt || "";
    var copyLabel = config.copyLabel || "Copy for AI";
    var canonicalBase = config.canonicalBase || "";
    var pagePath = config.pagePath || "";

    // Resolve the full markdown URL for AI provider links
    // Use current page location to get the correct deployed path (includes dev/, stable/, etc.)
    var fullMdUrl;
    if (mdUrl) {
      var loc = window.location;
      var basePath = loc.pathname.replace(/[^/]*$/, "");
      fullMdUrl = loc.origin + basePath + mdUrl.replace(/^\.\//, "");
    } else {
      fullMdUrl = "";
    }

    // --- Build DOM ---

    var container = document.createElement("div");
    container.className = "copybutton-container";

    // Dropdown (positioned above the button group)
    var dropdown = document.createElement("div");
    dropdown.className = "copybutton-dropdown dropdown-content";
    dropdown.setAttribute("role", "menu");
    container.appendChild(dropdown);

    // Tooltip
    var tooltip = document.createElement("span");
    tooltip.className = "copybutton-tooltip";
    tooltip.textContent = "Copied!";
    container.appendChild(tooltip);

    // Button group
    var group = document.createElement("div");
    group.className = "copybutton-group";

    // --- Copy for AI button ---
    var copyBtn = document.createElement("button");
    copyBtn.className = "button";
    copyBtn.title = "Copy page as Markdown for AI";
    copyBtn.setAttribute("aria-label", "Copy page as Markdown for AI");
    copyBtn.innerHTML =
      '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">' +
      '<rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>' +
      '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>' +
      '</svg>' +
      '<span>' + copyLabel + '</span>';

    copyBtn.addEventListener("click", function () {
      fetchAndCopy(mdUrl, function (err) {
        if (err) {
          console.error("DocumenterCopyButton: copy failed", err);
          return;
        }
        showTooltip(tooltip);
      });
    });

    group.appendChild(copyBtn);

    // --- Open in... button ---
    if (providers.length > 0) {
      var openBtn = document.createElement("button");
      openBtn.className = "button";
      openBtn.title = "Open in AI";
      openBtn.setAttribute("aria-label", "Open in AI");
      openBtn.innerHTML =
        '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">' +
        '<circle cx="12" cy="12" r="10"></circle>' +
        '<line x1="2" y1="12" x2="22" y2="12"></line>' +
        '<path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"></path>' +
        '</svg>' +
        '<span>AI</span>' +
        '<svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round">' +
        '<polyline points="6 9 12 15 18 9"></polyline>' +
        '</svg>';

      // Build dropdown items
      providers.forEach(function (provider) {
        var item = document.createElement("a");
        item.className = "dropdown-item";
        item.setAttribute("role", "menuitem");
        var nameSpan = document.createElement("span");
        nameSpan.textContent = provider.name || "Open";
        item.appendChild(nameSpan);
        var redirectSvg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        redirectSvg.setAttribute("width", "8");
        redirectSvg.setAttribute("height", "8");
        redirectSvg.setAttribute("viewBox", "0 0 8 8");
        redirectSvg.setAttribute("fill", "none");
        redirectSvg.setAttribute("aria-hidden", "true");
        redirectSvg.classList.add("copybutton-redirect-icon");
        var redirectPath = document.createElementNS("http://www.w3.org/2000/svg", "path");
        redirectPath.setAttribute("d", "M2.91667 0.75H7.25M7.25 0.75V5.08333M7.25 0.75L0.75 7.25");
        redirectPath.setAttribute("stroke", "currentColor");
        redirectPath.setAttribute("stroke-width", "1.5");
        redirectPath.setAttribute("stroke-linecap", "round");
        redirectPath.setAttribute("stroke-linejoin", "round");
        redirectSvg.appendChild(redirectPath);
        item.appendChild(redirectSvg);
        item.href = "#";
        item.addEventListener("click", function (e) {
          e.preventDefault();
          e.stopPropagation();
          var query = prompt ? prompt + " " + fullMdUrl : fullMdUrl;
          var url =
            (provider.base || "") +
            encodeURIComponent(query) +
            (provider.suffix || "");
          window.open(url, "_blank", "noopener,noreferrer");
          closeDropdown();
        });
        dropdown.appendChild(item);
      });

      openBtn.addEventListener("click", function (e) {
        e.stopPropagation();
        var isOpen = dropdown.classList.contains("is-open");
        closeDropdown();
        if (!isOpen) {
          openBtn.classList.add("active");
          dropdown.classList.add("is-open");
        }
      });

      group.appendChild(openBtn);
    }

    container.appendChild(group);
    document.body.appendChild(container);

    // Dismiss dropdown on outside click
    document.addEventListener("click", function () {
      closeDropdown();
    });

    // --- Helpers ---

    function fetchAndCopy(url, callback) {
      // Try embedded content first (works with file:// and offline)
      var embedded = document.getElementById("documenter-copybutton-content");
      if (embedded) {
        writeToClipboard(embedded.textContent, callback);
        return;
      }
      // Fall back to fetching the .md file
      if (!url) {
        callback(new Error("No markdown URL configured"));
        return;
      }
      fetch(url)
        .then(function (response) {
          if (!response.ok) throw new Error("HTTP " + response.status);
          return response.text();
        })
        .then(function (text) {
          writeToClipboard(text, callback);
        })
        .catch(function (err) {
          callback(err);
        });
    }

    function writeToClipboard(text, callback) {
      if (navigator.clipboard && navigator.clipboard.writeText) {
        navigator.clipboard.writeText(text).then(
          function () { callback(null); },
          function () { fallbackCopy(text, callback); }
        );
      } else {
        fallbackCopy(text, callback);
      }
    }

    function fallbackCopy(text, callback) {
      var textarea = document.createElement("textarea");
      textarea.value = text;
      textarea.style.position = "fixed";
      textarea.style.top = "-9999px";
      textarea.style.left = "-9999px";
      document.body.appendChild(textarea);
      textarea.focus();
      textarea.select();
      var success = false;
      try {
        success = document.execCommand("copy");
      } catch (e) {
        // ignore
      }
      document.body.removeChild(textarea);
      success ? callback(null) : callback(new Error("execCommand copy failed"));
    }

    function showTooltip(el) {
      el.classList.add("is-visible");
      setTimeout(function () {
        el.classList.remove("is-visible");
      }, 1800);
    }

    function closeDropdown() {
      dropdown.classList.remove("is-open");
      var active = container.querySelectorAll(".button.active");
      active.forEach(function (btn) { btn.classList.remove("active"); });
    }
  });
})();
