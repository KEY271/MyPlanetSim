import { createServer } from "node:http";

const server = createServer((request, response) => {
  response.setHeader("content-type", "application/json; charset=utf-8");
  response.statusCode = 404;
  response.end(JSON.stringify({ error: "gateway endpoints are not enabled yet" }));
});

server.listen(0, "127.0.0.1", () => {
  const address = server.address();
  console.log(`gateway listening on ${typeof address === "object" ? address?.port : address}`);
});
